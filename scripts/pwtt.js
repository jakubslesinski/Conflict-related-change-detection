/*
==================================================================================
📌 PAIRWISE WELCH'S T-TEST (PWTT) for SAR Change Detection
==================================================================================
Description:
- This script implements the Pairwise Welch's T-Test for detecting changes in 
  Sentinel-1 SAR time series data
- PWTT compares pre-event and post-event SAR backscatter statistics to identify
  statistically significant changes, particularly useful for conflict monitoring
- The algorithm applies Lee filtering for speckle reduction and uses multiple
  urban masks to focus analysis on built-up areas
- Outputs include T-statistic maps and binary damage/change classifications

Based on: https://github.com/oballinger/PWTT
Author: Based on Oliver Ballinger's implementation
References:
- Ballinger, O. et al. (2023). "Monitoring war destruction using high-resolution 
  satellite imagery and deep learning"
- Welch, B. L. (1947). "The generalization of Student's problem when several 
  different population variances are involved"
- Lee, J. S. (1980). "Digital image enhancement and noise filtering by use of 
  local statistics"
==================================================================================
*/

// =========================================================================
// PWTT (Pairwise Welch's T-Test) for Donbas Region - Modified Version
// Based on: https://github.com/oballinger/PWTT
// =========================================================================

// =========================================================================
// QUICK START OPTION
// =========================================================================
// If you're getting errors with MS Buildings or other data, set this to true:
var USE_SIMPLE_MODE = false;  // Set to true for simple mode with Dynamic World only
// =========================================================================

// =========================================================================
// REGION CONFIGURATION - DONBAS
// =========================================================================

// Assuming Donbas is already defined as a FeatureCollection
// If not, define it here:
// var Donbas = ee.FeatureCollection('path/to/donbas/boundaries');

// Check region size and complexity
var regionArea = Donbas.geometry().area();
var regionAreaKm2 = regionArea.divide(1000000); // Convert to km²
print('📐 Donbas region area (km²):', regionAreaKm2);

// If region is too large, consider simplifying
var maxAreaKm2 = 50000; // Maximum recommended area
var isRegionTooLarge = regionAreaKm2.gt(maxAreaKm2);
print('⚠️ Region might be too large for processing:', isRegionTooLarge);

// Simplify geometry if needed to reduce processing complexity
var simplifiedDonbas = Donbas.map(function(feature) {
  return feature.simplify(100); // Simplify to 100m tolerance
});

// Study area: Donbas region
var geometry = simplifiedDonbas.geometry();

// Set default basemap to satellite
Map.setOptions('satellite');

// Center the map on the Donbas region
Map.centerObject(simplifiedDonbas, 10);
Map.addLayer(simplifiedDonbas.style({color: 'red', fillColor: '00000000'}), {}, 'Donbas Region');

// =========================================================================
// PWTT PARAMETERS
// =========================================================================

// War/conflict start dates
var warStart = '2022-01-01';  // Start of conflict
var inferenceStart = '2025-01-01';  // Analysis start date
var preInterval = 12;  // Months before war start for baseline (inferenceStart-12 months)
var postInterval = 6;  // Months after inference start for analysis
//var warStart = '2022-02-24';
//var inferenceStart = '2022-07-01';  // Start of second half of year
//var preInterval = 12;               // Year before invasion
//var postInterval = 6;               // July-December 2022

// Thresholding parameters
var enableThresholding = true;
var tStatThreshold = 3;  // T-statistic threshold for change detection
var damageThreshold = 3;  // Damage classification threshold
var urbanThreshold = 0.1;  // Urban area threshold from Dynamic World

// Export parameters
var exportScale = 10;
var exportDir = 'PWTT_Donbas_Export';
var gridScale = 500;

// Orbit mode
var orbitMode = 'ASCENDING';

// Active mask selection (change this to switch between masks)
var activeMask = 'DYNAMIC_WORLD'; // Options: 'DYNAMIC_WORLD', 'MS_BUILDINGS', 'ESA_WORLDCOVER', 'ALL'

// =========================================================================
// LEE FILTER IMPLEMENTATION
// =========================================================================

var leeFilter = function(image) {
  var KERNEL_SIZE = 2;
  var bandNames = image.bandNames().remove('angle');
  
  // S1-GRD images are multilooked 5 times in range
  var enl = 5;
  
  // Compute the speckle standard deviation
  var eta = 1.0 / Math.sqrt(enl);
  eta = ee.Image.constant(eta);
  
  // MMSE estimator
  var oneImg = ee.Image.constant(1);
  
  var reducers = ee.Reducer.mean().combine({
    reducer2: ee.Reducer.variance(),
    sharedInputs: true
  });
  
  var stats = image.select(bandNames).reduceNeighborhood({
    reducer: reducers,
    kernel: ee.Kernel.square(KERNEL_SIZE / 2, 'pixels'),
    optimization: 'window'
  });
  
  var meanBand = bandNames.map(function(bandName) {
    return ee.String(bandName).cat('_mean');
  });
  var varBand = bandNames.map(function(bandName) {
    return ee.String(bandName).cat('_variance');
  });
  
  var zBar = stats.select(meanBand);
  var varz = stats.select(varBand);
  
  // Estimate weight
  var varx = varz.subtract(zBar.pow(2).multiply(eta.pow(2)))
    .divide(oneImg.add(eta.pow(2)));
  var b = varx.divide(varz);
  
  // If b is negative, set it to zero
  var newB = b.where(b.lt(0), 0);
  
  var output = oneImg.subtract(newB).multiply(zBar.abs())
    .add(newB.multiply(image.select(bandNames)));
  output = output.rename(bandNames);
  
  return image.addBands(output, null, true);
};

// =========================================================================
// T-TEST IMPLEMENTATION
// =========================================================================

var ttest = function(s1, inferenceStartDate, warStartDate, preInt, postInt) {
  // Convert dates to ee.Date objects
  var infStart = ee.Date(inferenceStartDate);
  var warStart = ee.Date(warStartDate);
  
  // Filter pre-event period
  var pre = s1.filterDate(
    warStart.advance(ee.Number(preInt).multiply(-1), "month"),
    warStart
  );
  
  // Filter post-event period  
  var post = s1.filterDate(
    infStart, 
    infStart.advance(postInt, "month")
  );
  
  // Calculate statistics for pre-event period
  var preMean = pre.mean();
  var preSd = pre.reduce(ee.Reducer.stdDev());
  var preN = ee.Number(pre.aggregate_array('orbitNumber_start')
    .distinct().size());
  
  // Calculate statistics for post-event period
  var postMean = post.mean();
  var postSd = post.reduce(ee.Reducer.stdDev());
  var postN = ee.Number(post.aggregate_array('orbitNumber_start')
    .distinct().size());
  
  // Calculate pooled standard deviation
  var pooledSd = preSd.pow(2)
    .multiply(preN.subtract(1))
    .add(postSd.pow(2).multiply(postN.subtract(1)))
    .divide(preN.add(postN).subtract(2))
    .sqrt();
  
  // Calculate denominator of t-test
  var denom = pooledSd.multiply(
    ee.Image(1).divide(preN)
      .add(ee.Image(1).divide(postN))
      .sqrt()
  );
  
  // Calculate t-test
  var change = postMean.subtract(preMean)
    .divide(denom).abs();
  
  return change;
};

// =========================================================================
// URBAN MASKS CREATION
// =========================================================================

print('🏗️ Creating urban masks from different sources...');

// Convert dates
var inferenceStartDate = ee.Date(inferenceStart);
var warStartDate = ee.Date(warStart);

// 1. Google Dynamic World mask
var dynamicWorldMask = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1')
  .filterDate(warStartDate.advance(-1 * preInterval, 'months'), warStartDate)
  .select('built')
  .mean()
  .gt(urbanThreshold)
  .rename('dynamic_world_mask')
  .clip(simplifiedDonbas);

// 2. Microsoft Building Footprints mask - with better error handling
var msBuildingsMask;
var msAvailable = false;
try {
  var msBuildings = ee.FeatureCollection('projects/sat-io/open-datasets/MSBuildings/Ukraine');
  
  // Check if collection exists and has features
  var msCount = msBuildings.size();
  
  msBuildings = msBuildings.filterBounds(simplifiedDonbas);
  
  // Convert footprints to raster
  msBuildingsMask = msBuildings
    .map(function(feat) { return feat.set('presence', 1); })
    .reduceToImage(['presence'], ee.Reducer.first())
    .gt(0)
    .rename('ms_buildings_mask')
    .reproject({crs: 'EPSG:4326', scale: 10})
    .clip(simplifiedDonbas);
    
  msAvailable = true;
  print('✅ Microsoft Building Footprints loaded successfully');
  print('📊 MS Buildings in region:', msBuildings.size());
} catch (e) {
  print('⚠️ Microsoft Building Footprints not available - using fallback');
  print('💡 Error:', e);
  msBuildingsMask = ee.Image(0).rename('ms_buildings_mask').clip(simplifiedDonbas);
  msAvailable = false;
}

// 3. ESA WorldCover mask
var esaWorldCoverMask = ee.ImageCollection("ESA/WorldCover/v200")
  .first()
  .select('Map')
  .eq(50)  // 50 = built-up areas
  .rename('esa_worldcover_mask')
  .reproject({crs: 'EPSG:4326', scale: 10})
  .clip(simplifiedDonbas);

// Combine available masks for comparison
var allMasks = ee.Image.cat([
  dynamicWorldMask,
  msBuildingsMask,
  esaWorldCoverMask
]);

// =========================================================================
// MAIN PWTT PROCESSING
// =========================================================================

print('🚀 Starting PWTT processing for Donbas region...');

// Get unique orbit numbers
var orbits = ee.ImageCollection("COPERNICUS/S1_GRD_FLOAT")
  .filter(ee.Filter.listContains("transmitterReceiverPolarisation", "VH"))
  .filter(ee.Filter.eq("instrumentMode", "IW"))
  .filter(ee.Filter.eq("orbitProperties_pass", orbitMode))
  .filterBounds(simplifiedDonbas)
  .filter(ee.Filter.contains('.geo', simplifiedDonbas.geometry()))
  .filterDate(inferenceStartDate, inferenceStartDate.advance(postInterval, 'months'))
  .aggregate_array('relativeOrbitNumber_start')
  .distinct();

print('📡 Number of unique orbits:', orbits.size());
print('📡 Orbit numbers:', orbits);

// Count the total number of Sentinel-1 images used in analysis
var sentinelImages = ee.ImageCollection("COPERNICUS/S1_GRD_FLOAT")
  .filter(ee.Filter.listContains("transmitterReceiverPolarisation", "VH"))
  .filter(ee.Filter.eq("instrumentMode", "IW"))
  .filter(ee.Filter.eq("orbitProperties_pass", orbitMode))
  .filterBounds(simplifiedDonbas)
  .filter(ee.Filter.contains('.geo', simplifiedDonbas.geometry()))
  .filterDate(warStartDate.advance(-preInterval, 'months'),
              inferenceStartDate.advance(postInterval, 'months'));

print('🛰️ Number of Sentinel-1 images used in analysis:', sentinelImages.size());

// Function to process image with specific mask
var processWithMask = function(image, mask, maskName) {
  var masked = image.updateMask(mask);
  
  // Add convolution layers
  var k50 = masked.convolve(ee.Kernel.circle(50, 'meters', true)).rename('k50');
  var k100 = masked.convolve(ee.Kernel.circle(100, 'meters', true)).rename('k100');
  var k150 = masked.convolve(ee.Kernel.circle(150, 'meters', true)).rename('k150');
  
  // Add damage band
  var damage = masked.select('max_change').gt(damageThreshold).rename('damage');
  var result = masked.addBands(damage);
  result = result.addBands([k50, k100, k150]);
  
  // Calculate T_statistic
  result = result.addBands(
    masked.select('max_change')
      .add(k50)
      .add(k100)
      .add(k150)
      .divide(4)
      .rename('T_statistic')
  );
  
  return result.select('T_statistic', 'damage').toFloat()
    .set('mask_type', maskName);
};

// =========================================================================
// QUICK START - SIMPLIFIED VERSION WITH DYNAMIC WORLD ONLY
// =========================================================================

if (USE_SIMPLE_MODE) {
  print('🚀 Running in SIMPLE MODE - using only Dynamic World mask');
  
  // Map over each orbit - fixed to avoid client-side operations
  var mapOrbit = function(orbit) {
    var s1 = ee.ImageCollection("COPERNICUS/S1_GRD_FLOAT")
      .filter(ee.Filter.listContains("transmitterReceiverPolarisation", "VH"))
      .filter(ee.Filter.eq("instrumentMode", "IW"))
      .filter(ee.Filter.eq("relativeOrbitNumber_start", orbit))
      .filter(ee.Filter.eq("orbitProperties_pass", orbitMode))
      .map(leeFilter)
      .select(['VV', 'VH'])
      .map(function(image) {
        return image.log();
      })
      .filterBounds(simplifiedDonbas);
    
    var ttestResult = ttest(s1, inferenceStartDate, warStartDate, preInterval, postInterval);
    return ttestResult;
  };
  
  // Process all orbits
  var ttestImages = ee.ImageCollection(orbits.map(mapOrbit));
  print('🔍 Number of t-test results:', ttestImages.size());
  
  // Create base image - max across all orbits
  var baseImage = ttestImages.max();
  
  // Add max change band
  baseImage = baseImage.addBands(
    baseImage.select('VV').max(baseImage.select('VH')).rename('max_change')
  ).select('max_change');
  
  // Apply focal median
  baseImage = baseImage.focalMedian(10, 'gaussian', 'meters')
    .clip(simplifiedDonbas);
  
  // Process only with Dynamic World
  var result = processWithMask(baseImage, dynamicWorldMask, 'Dynamic World');
  
  // Color palette for T-statistic
  var tStatPalette = ['yellow', 'orange', 'red', 'purple'];
  
  // Visualization
  Map.addLayer(dynamicWorldMask.selfMask(), {palette: ['cyan']}, 'Urban Mask: Dynamic World', false);
  
  Map.addLayer(result.select('T_statistic'), {
    min: 3, max: 5, palette: tStatPalette
  }, 'PWTT Result (Dynamic World)', true);
  
  Map.addLayer(result.select('damage').selfMask(), {
    palette: ['red']
  }, 'Damage Areas', false);
  
  // Simple statistics
  var changePixels = result.select('damage').reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: simplifiedDonbas.geometry(),
    scale: 100,
    maxPixels: 1e9
  });
  
  print('📊 Damage statistics:', changePixels);
  
  // Export
  Export.image.toDrive({
    image: result.select('T_statistic').multiply(1000).toInt16(),
    description: 'PWTT_TStatistic_DynamicWorld_Donbas_Simple',
    folder: exportDir,
    region: simplifiedDonbas.geometry(),
    scale: exportScale,
    crs: 'EPSG:4326',
    maxPixels: 1e13
  });
  
  print('✅ Simple mode processing complete! Check the Tasks tab for export.');
  
} else {
  // Continue with full multi-mask processing...

// Map over each orbit - fixed to avoid client-side operations
var mapOrbit = function(orbit) {
  var s1 = ee.ImageCollection("COPERNICUS/S1_GRD_FLOAT")
    .filter(ee.Filter.listContains("transmitterReceiverPolarisation", "VH"))
    .filter(ee.Filter.eq("instrumentMode", "IW"))
    .filter(ee.Filter.eq("relativeOrbitNumber_start", orbit))
    .filter(ee.Filter.eq("orbitProperties_pass", orbitMode))
    .map(leeFilter)
    .select(['VV', 'VH'])
    .map(function(image) {
      return image.log();
    })
    .filterBounds(simplifiedDonbas);
  
  var ttestResult = ttest(s1, inferenceStartDate, warStartDate, preInterval, postInterval);
  return ttestResult;
};

// Process all orbits
var ttestImages = ee.ImageCollection(orbits.map(mapOrbit));
print('🔍 Number of t-test results:', ttestImages.size());

// Create base image - max across all orbits
var baseImage = ttestImages.max();

// Add max change band
baseImage = baseImage.addBands(
  baseImage.select('VV').max(baseImage.select('VH')).rename('max_change')
).select('max_change');

// Apply focal median
baseImage = baseImage.focalMedian(10, 'gaussian', 'meters')
  .clip(simplifiedDonbas);

// Process with each mask
var dynamicWorldResult = processWithMask(baseImage, dynamicWorldMask, 'Dynamic World');
var msBuildingsResult = processWithMask(baseImage, msBuildingsMask, 'MS Buildings');
var esaWorldCoverResult = processWithMask(baseImage, esaWorldCoverMask, 'ESA WorldCover');

// =========================================================================
// VISUALIZATION
// =========================================================================

// Color palette for T-statistic
var tStatPalette = ['yellow', 'orange', 'red', 'purple'];
var maskPalette = ['black', 'white'];

// Add urban mask layers
Map.addLayer(dynamicWorldMask.selfMask(), {palette: ['cyan']}, 'Mask: Dynamic World', false);
Map.addLayer(msBuildingsMask.selfMask(), {palette: ['yellow']}, 'Mask: MS Buildings', false);
Map.addLayer(esaWorldCoverMask.selfMask(), {palette: ['magenta']}, 'Mask: ESA WorldCover', false);

// Add PWTT results for each mask
Map.addLayer(dynamicWorldResult.select('T_statistic'), {
  min: 3, max: 5, palette: tStatPalette
}, 'PWTT: Dynamic World', true);

Map.addLayer(msBuildingsResult.select('T_statistic'), {
  min: 3, max: 5, palette: tStatPalette
}, 'PWTT: MS Buildings', false);

Map.addLayer(esaWorldCoverResult.select('T_statistic'), {
  min: 3, max: 5, palette: tStatPalette
}, 'PWTT: ESA WorldCover', false);

// =========================================================================
// MASK COMPARISON STATISTICS
// =========================================================================

print('📊 Computing mask statistics...');

// Function to calculate mask statistics - fixed
var calculateMaskStats = function(mask, name) {
  var bandName = mask.bandNames().get(0);
  
  var stats = mask.reduceRegion({
    reducer: ee.Reducer.sum().combine(ee.Reducer.count(), '', true),
    geometry: simplifiedDonbas.geometry(),
    scale: 100,
    maxPixels: 1e9
  });
  
  // Get the keys from the stats dictionary
  var sumKey = ee.String(bandName).cat('_sum');
  var countKey = ee.String(bandName).cat('_count');
  
  var maskArea = ee.Number(stats.get(sumKey))
    .multiply(100 * 100 / 1000000); // Convert to km²
  
  var totalPixels = stats.get(countKey);
  var maskPercentage = ee.Number(stats.get(sumKey))
    .divide(totalPixels).multiply(100);
  
  return ee.Feature(null, {
    'mask_name': name,
    'urban_area_km2': maskArea,
    'urban_percentage': maskPercentage
  });
};

// Calculate statistics for each mask
var maskStats = ee.FeatureCollection([
  calculateMaskStats(dynamicWorldMask, 'Dynamic World'),
  calculateMaskStats(msBuildingsMask, 'MS Buildings'),
  calculateMaskStats(esaWorldCoverMask, 'ESA WorldCover')
]);

print('📊 Urban Mask Statistics:', maskStats);

// =========================================================================
// INTERACTIVE COMPARISON
// =========================================================================

// Create a panel for interactive comparison
Map.style().set('cursor', 'crosshair');

var panel = ui.Panel();
panel.style().set({
  width: '500px',
  position: 'bottom-left'
});

Map.add(panel);

// Register click handler
Map.onClick(function(coords) {
  panel.clear();
  var point = ee.Geometry.Point(coords.lon, coords.lat);
  
  // Check if point is within Donbas
  var withinDonbas = simplifiedDonbas.geometry().contains(point);
  var isInside = withinDonbas.getInfo();
  
  if (!isInside) {
    panel.add(ui.Label('⚠️ Click within the Donbas region to see comparisons'));
    return;
  }
  
  panel.add(ui.Label('📍 Location: ' + coords.lon.toFixed(4) + ', ' + coords.lat.toFixed(4)));
  panel.add(ui.Label(''));
  
  // Get mask values at point
  var maskValues = allMasks.reduceRegion({
    reducer: ee.Reducer.first(),
    geometry: point,
    scale: 10
  });
  
  panel.add(ui.Label('🏗️ Urban Mask Values:', {fontWeight: 'bold'}));
  panel.add(ui.Label('Dynamic World: ' + maskValues.get('dynamic_world_mask').getInfo()));
  panel.add(ui.Label('MS Buildings: ' + maskValues.get('ms_buildings_mask').getInfo()));
  panel.add(ui.Label('ESA WorldCover: ' + maskValues.get('esa_worldcover_mask').getInfo()));
  panel.add(ui.Label(''));
  
  // Get T-statistic values for each mask variant
  var getTStat = function(image, name) {
    var value = image.select('T_statistic').reduceRegion({
      reducer: ee.Reducer.first(),
      geometry: point,
      scale: 10
    }).get('T_statistic');
    
    return name + ': ' + ee.Number(value).format('%.3f').getInfo();
  };
  
  panel.add(ui.Label('📊 T-statistic by Mask:', {fontWeight: 'bold'}));
  panel.add(ui.Label(getTStat(dynamicWorldResult, 'Dynamic World')));
  panel.add(ui.Label(getTStat(msBuildingsResult, 'MS Buildings')));
  panel.add(ui.Label(getTStat(esaWorldCoverResult, 'ESA WorldCover')));
  
  // Create time series chart
  var s1Point = ee.ImageCollection("COPERNICUS/S1_GRD_FLOAT")
    .filter(ee.Filter.listContains("transmitterReceiverPolarisation", "VH"))
    .filter(ee.Filter.eq("instrumentMode", "IW"))
    .filter(ee.Filter.eq("orbitProperties_pass", orbitMode))
    .filterBounds(point)
    .select(['VV', 'VH'])
    .filterDate(warStartDate.advance(-preInterval, 'months'), 
                inferenceStartDate.advance(postInterval, 'months'));
  
  var chart = ui.Chart.image.series({
    imageCollection: s1Point,
    region: point,
    reducer: ee.Reducer.mean(),
    scale: 10
  }).setOptions({
    title: 'Sentinel-1 Time Series at Point',
    hAxis: {title: 'Date'},
    vAxis: {title: 'Backscatter (dB)'},
    lineWidth: 2,
    pointSize: 4,
    series: {
      0: {color: 'blue'},  // VV
      1: {color: 'red'}    // VH
    }
  });
  
  panel.add(chart);
});

// =========================================================================
// EXPORTS
// =========================================================================

var dateString = warStartDate.format('YYYY-MM-dd').getInfo() + '_to_' + 
                 inferenceStartDate.format('YYYY-MM-dd').getInfo();

print('📤 Creating export tasks for mask variants...');

// Function to export results for a specific mask
var exportMaskResults = function(result, maskName) {
  var cleanName = maskName.replace(/\s+/g, '_').replace(/[()]/g, '');
  
  // Export T-statistic
  Export.image.toDrive({
    image: result.select('T_statistic').multiply(1000).toInt16(),
    description: 'PWTT_TStatistic_' + cleanName + '_Donbas_' + dateString,
    folder: exportDir,
    fileNamePrefix: 'PWTT_TStatistic_' + cleanName + '_Donbas_' + dateString,
    region: simplifiedDonbas.geometry(),
    scale: exportScale,
    crs: 'EPSG:4326',
    maxPixels: 1e13
  });
  
  // Export damage mask
  Export.image.toDrive({
    image: result.select('damage').toByte(),
    description: 'PWTT_Damage_' + cleanName + '_Donbas_' + dateString,
    folder: exportDir,
    fileNamePrefix: 'PWTT_Damage_' + cleanName + '_Donbas_' + dateString,
    region: simplifiedDonbas.geometry(),
    scale: exportScale,
    crs: 'EPSG:4326',
    maxPixels: 1e13
  });
};

// Export results for each mask variant
exportMaskResults(dynamicWorldResult, 'Dynamic_World');
exportMaskResults(msBuildingsResult, 'MS_Buildings');
exportMaskResults(esaWorldCoverResult, 'ESA_WorldCover');

// Export mask comparison statistics
Export.table.toDrive({
  collection: maskStats,
  description: 'PWTT_MaskStatistics_Donbas_' + dateString,
  folder: exportDir,
  fileFormat: 'CSV'
});

// Export all masks as multi-band image
Export.image.toDrive({
  image: allMasks.toByte(),
  description: 'PWTT_AllMasks_Donbas_' + dateString,
  folder: exportDir,
  fileNamePrefix: 'PWTT_AllMasks_Donbas_' + dateString,
  region: simplifiedDonbas.geometry(),
  scale: exportScale,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

print('📤 Export tasks created. Check the Tasks tab to run exports.');

// =========================================================================
// SUMMARY
// =========================================================================

var summary = ee.Dictionary({
  '🌍 Region': 'Donbas',
  '📅 War Start': warStart,
  '📅 Analysis Date': inferenceStart,
  '📊 Pre-interval': preInterval + ' months',
  '📊 Post-interval': postInterval + ' months',
  '📡 Orbit Mode': orbitMode,
  '🎚️ T-stat Threshold': tStatThreshold,
  '🏚️ Damage Threshold': damageThreshold,
  '🏙️ Urban Threshold': urbanThreshold,
  '📏 Export Scale': exportScale + ' meters',
  '🎭 Mask Variants': '3 (Dynamic World, MS Buildings, ESA WorldCover)'
});
print('📊 PWTT Analysis Summary:', summary);

// =========================================================================
// USAGE INSTRUCTIONS
// =========================================================================

print('');
print('🎛️ PWTT PARAMETERS:');
print('• 3 urban mask variants are processed');
print('• Click on map to compare T-statistics across masks');
print('• Each mask produces separate export files');
print('');
print('🎭 MASK DESCRIPTIONS:');
print('• Dynamic World: Google/Sentinel-2 based ML classification');
print('• MS Buildings: Microsoft building footprints (vectorized)');
print('• ESA WorldCover: ESA 10m land cover classification');
print('');
print('📊 INTERPRETATION:');
print('• Different masks may show different damage patterns');
print('• Compare results to choose optimal mask for your analysis');
print('• Higher T-statistics indicate stronger statistical significance');
print('');
print('🔬 STATISTICAL METHODOLOGY:');
print('• Welch\'s t-test compares pre/post-event SAR statistics');
print('• Lee filtering reduces speckle noise in SAR data');
print('• Spatial convolution provides contextual smoothing');
print('• Urban masks focus analysis on built-up areas');
print('• Log transformation normalizes SAR backscatter values');
print('');
print('💡 TIPS:');
print('• Toggle layers on/off to compare visually');
print('• Check mask statistics in console for coverage info');
print('• Export all variants for post-processing comparison');
print('• Use significance thresholds to control false positive rates');
print('');
print('✅ PWTT processing complete!');
} // End of full processing mode

// =========================================================================
// ALGORITHM PERFORMANCE NOTES
// =========================================================================
print('');
print('🔧 PERFORMANCE CONSIDERATIONS:');
print('• Processing time depends on region size and time series length');
print('• Lee filter adds computational overhead but improves results');
print('• Spatial convolution kernels provide noise reduction');
print('• Multiple masks allow sensitivity analysis of results');
print('• Export scale affects file size and processing time');
print('');
print('📚 ALGORITHM REFERENCES:');
print('• PWTT methodology: Oliver Ballinger et al.');
print('• Welch\'s t-test: B.L. Welch (1947)');
print('• Lee speckle filter: J.S. Lee (1980)');
print('• Implementation: Google Earth Engine platform');