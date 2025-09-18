/*
==================================================================================
📌 FROZEN BACKGROUND DETECTION using a Sliding Window on Sentinel-1 Time Series
==================================================================================
Description:
- This script calculates a frozen background on a Sentinel-1 time series  
- A frozen background is an average image of a scene in which 
anomalies identified using the coefficient of variation have been removed
- The script processes the two polarisations (VV and VH) separately.
- In this implementation, instead of using an iterative approach, it directly 
operates on the sorted list of pixel values to calculate cumulative coefficients
of variation (CV) and determines the number of values that remain below the threshold
Author: Elise Colin (Onera) from:
Reference:
Taillade, T.; Thirion-Lefevre, L.; Guinvarc'h, R. 
Detecting Ephemeral Objects in SAR Time-Series Using Frozen Background-Based 
Change Detection. Remote Sens. 2020, 12, 1720. https://doi.org/10.3390/rs12111720 

==================================================================================
*/

// =========================================================================
// 1. User parameters 
// =========================================================================
var str2='2022-01-01';     // End of the Observation 'YY-MM-dd'
var durationMonths = 12;    // Duration of the time series ( in month)

var mu_threshold = 0.2286;

// Check region size and complexity
var regionArea = Donbas.geometry().area();
var regionAreaKm2 = regionArea.divide(1000000); // Convert to km²
print('📐 Donbas region area (km²):', regionAreaKm2);

// If region is too large, consider simplifying or processing in tiles
var maxAreaKm2 = 50000; // Maximum recommended area
var isRegionTooLarge = regionAreaKm2.gt(maxAreaKm2);
print('⚠️ Region might be too large for processing:', isRegionTooLarge);

// Simplify geometry if needed to reduce processing complexity
var simplifiedDonbas = Donbas.map(function(feature) {
  return feature.simplify(100); // Simplify to 100m tolerance
});

// -------------------------------------------------------------------------
// COMPUTE COLLECTION
// -------------------------------------------------------------------------
var endDate = ee.Date(str2);
var startDate = endDate.advance(-durationMonths, 'month'); 

// Study area: Donbas region (use simplified version if original is too complex)
var geometry = simplifiedDonbas.geometry();

// Center the map on the Donbas region
Map.centerObject(simplifiedDonbas, 13);
//Map.addLayer(simplifiedDonbas, {color: 'red'}, 'Donbas Region', false);
Map.addLayer(simplifiedDonbas.style({color: 'red', fillColor: '00000000'}));

// =========================================================================
// 2. Loading and initial filtering of the Sentinel-1 collection
// =========================================================================
var s1Collection = ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT')
  .filterBounds(simplifiedDonbas);  // Filter by Donbas region

var NbOrbit = s1Collection.aggregate_count_distinct('relativeOrbitNumber_start');
var ListOrbits = s1Collection.aggregate_array('relativeOrbitNumber_start');
// find orbit numbers and their frequency
var freq = ee.Dictionary(ee.List(ListOrbits).reduce(ee.Reducer.frequencyHistogram()));
var array = ee.Array([freq.keys().map(ee.Number.parse), freq.values()]);
// orbit choice : first, the one with the max frequency
var frequences = array.slice(0,-1);
var arraysort = array.sort(frequences);
var index = ee.Number(NbOrbit).add(-1);
var orbite = arraysort.get([0,ee.Number(index)]);
// find images with the chosen orbit
s1Collection = s1Collection.filterMetadata('relativeOrbitNumber_start', 'equals', orbite);

// Sort Sentinel-1 images by ascending date
var sortedCollection = s1Collection.sort('system:time_start', true);
// Find the last date before endDate
var beforeEnd = sortedCollection.filter(ee.Filter.date(startDate, endDate))
                                .sort('system:time_start', false)
                                .first()
                                .get('system:time_start');
var realEndDate = ee.Date(beforeEnd);
// Find the first date after realEndDate
var afterEnd = sortedCollection.filter(ee.Filter.date(endDate, endDate.advance(1, 'month')))
                               .sort('system:time_start', true)
                               .first()
                               .get('system:time_start');
var shiftedEndDate = ee.Date(afterEnd);

var analysisSummary = ee.Dictionary({
  '📅 Target End Date': endDate.format('YYYY-MM-dd'),
  '✅ Actual Last Available Date': realEndDate.format('YYYY-MM-dd'),
  '📷 Date of Analyzed Image': shiftedEndDate.format('YYYY-MM-dd')
});
print('🛰️ Temporal Analysis Summary', analysisSummary);

// Filter Sentinel-1 collections with the actual dates found
var s1_1 = s1Collection.filterDate(startDate, realEndDate.advance(1, 'day'));  // Main series
var s1_2 = s1Collection.filterDate(shiftedEndDate, shiftedEndDate.advance(1, 'day'));  // Shifted series

var s1Collection = s1Collection.filterDate(startDate, shiftedEndDate.advance(1, 'day')); // Main series

// Print collection information for debugging
print('📊 S1_1 collection size:', s1_1.size());
print('📊 S1_2 collection size:', s1_2.size());
print('📊 Total collection size:', s1Collection.size());

// Check if we have enough images
var minImages = 5; // Minimum number of images needed for frozen background
var s1_1_size = s1_1.size();
var hasEnoughImages = s1_1_size.gte(minImages);

print('✅ Has enough images for analysis:', hasEnoughImages);

var Inew = s1_2.first();
var Inew_linearVV = Inew.select('VV').sqrt().rename('VV').clip(simplifiedDonbas);
var Inew_linearVH = Inew.select('VH').sqrt().rename('VH').clip(simplifiedDonbas);

// Convert Intensity float to amplitude: sqrt(intensity)
function convertToAmplitude(collection, polar) {
  return collection.select(polar).map(function(img) {
    return img.expression('sqrt(band)', { 'band': img.select(polar) })
              .rename(polar);
  });
}

// For comparison, compute the simple mean for each polarization
var meanVV = convertToAmplitude(s1_1, 'VV').mean().rename('meanVV');
var meanVH = convertToAmplitude(s1_1, 'VH').mean().rename('meanVH');

// =========================================================================
// 3. Frozen Background
// =========================================================================

/**
 * Simplified frozen background computation to avoid array slicing issues
 * Uses statistical approach based on temporal coefficient of variation
 */

// Alternative simpler method that avoids complex array operations
function computeFrozenBackgroundSimple(collection, threshold, polar) {
  var amplitudeCollection = convertToAmplitude(collection, polar);
  
  // Calculate temporal statistics
  var mean = amplitudeCollection.mean();
  var stdDev = amplitudeCollection.reduce(ee.Reducer.stdDev());
  var cv = stdDev.divide(mean.max(0.0001)); // Avoid division by zero
  
  // Use multiple percentiles for robustness
  var p10 = amplitudeCollection.reduce(ee.Reducer.percentile([10]));
  var p25 = amplitudeCollection.reduce(ee.Reducer.percentile([25]));
  var p50 = amplitudeCollection.reduce(ee.Reducer.percentile([50]));
  
  // Progressive selection based on CV threshold
  var frozenBg = mean
    .where(cv.gt(threshold), p50)
    .where(cv.gt(threshold * 1.5), p25)
    .where(cv.gt(threshold * 2), p10)
    .rename('frozenBackground_' + polar);
  
  return frozenBg;
}

/**
 * Original computeFrozenBackground function with fixes for array operations
 * Keeping it for reference but using the simpler version by default
 */
function computeFrozenBackground(collection, mu_threshold, polar) {
  var amplitudeCollection = convertToAmplitude(collection, polar); // Convert the collection to amplitude 
  
  // Check collection size and convert to integer
  var collectionSize = amplitudeCollection.size();
  var collectionSizeInt = collectionSize.toInt();
  print('Collection size for ' + polar + ':', collectionSize);
  
  // Step 1: Convert the collection to a per-pixel sorted array
  var pixelValuesArray = amplitudeCollection.toArray().arraySort();
  
  // Step 2: Compute the cumulative coefficient of variation (CV) per pixel
  var indices = ee.List.sequence(1, collectionSizeInt);
  var cumulativeCV = indices.map(function(i) {
    var iInt = ee.Number(i).toInt(); // Ensure integer type
    // Extract a sub-array with the first i elements - using integer indices
    var subArray = pixelValuesArray.arraySlice(0, 0, iInt);
    // Compute the mean and standard deviation of the sub-array
    var mean = subArray.arrayReduce(ee.Reducer.mean(), [0]);
    var stdDev = subArray.arrayReduce(ee.Reducer.stdDev(), [0]);
    
    // Handle division by zero
    var cv = stdDev.divide(mean.where(mean.eq(0), 1));
    
    // Rename the band using the current index
    var newName = ee.String('cv_').cat(ee.Number(i).toInt().format());
    return cv.rename(newName);
  });
  
  // Convert the list of CV images into an ImageCollection then into a multi-band image
  var cumulativeCVImage = ee.ImageCollection.fromImages(cumulativeCV).toBands();
  
  // Create a boolean image: true where CV is less than the threshold
  var belowThreshold = cumulativeCVImage.lt(mu_threshold);
  
  // Step 3: Count the number of iterations (bands) where CV is below the threshold
  var countBelow = belowThreshold.reduce(ee.Reducer.sum()).rename('countBelow').toInt();
  
  // Create an index array matching pixelValuesArray dimensions
  var sequentialIndices = ee.List.sequence(0, collectionSizeInt.subtract(1));
  var indexArray = ee.Image.constant(sequentialIndices).toArray(0).toInt();
  
  // Create mask where index < countBelow
  var mask = indexArray.lt(countBelow.toArray(0));
  
  // Apply mask and compute frozen background
  var maskedArray = pixelValuesArray.arrayMask(mask);
  var sumAmplitude = maskedArray.arrayReduce(ee.Reducer.sum(), [0]);
  var countNonZero = maskedArray.arrayMask(maskedArray).arrayReduce(ee.Reducer.count(), [0]);
  
  // Compute the frozen background as the mean of the retained amplitudes
  // Use maximum to avoid division by zero
  var frozenBackground = sumAmplitude.divide(countNonZero.max(1))
                                     .arrayProject([1])
                                     .arrayFlatten([['frozenBackground']])
                                     .rename('frozenBackground_' + polar);
  
  return frozenBackground;
}

// Compute the frozen background for both polarizations
// Using simplified approach to avoid array slicing issues
var frozenBackgroundVV = computeFrozenBackgroundSimple(s1_1, mu_threshold, 'VV');
var frozenBackgroundVH = computeFrozenBackgroundSimple(s1_1, mu_threshold, 'VH');

// Ensure proper band naming
frozenBackgroundVV = frozenBackgroundVV.select([0], ['frozenBackground_VV']);
frozenBackgroundVH = frozenBackgroundVH.select([0], ['frozenBackground_VH']);

// Clip all results to the Donbas region
frozenBackgroundVV = frozenBackgroundVV.clip(simplifiedDonbas);
frozenBackgroundVH = frozenBackgroundVH.clip(simplifiedDonbas);
meanVV = meanVV.clip(simplifiedDonbas);
meanVH = meanVH.clip(simplifiedDonbas);

// =========================================================================
// 4. Change Detection and Thresholding
// =========================================================================

// Threshold parameters
var changeThresholdVV = 0.15;  // Threshold for VV change (0-1 scale, adjust as needed)
var changeThresholdVH = 0.15;  // Threshold for VH change (0-1 scale, adjust as needed)
var minChangeArea = 10;        // Minimum connected pixels to consider as change
var changeType = 'both';       // 'increase', 'decrease', or 'both'

// Calculate difference between new image and frozen background
var diffVV = Inew_linearVV.subtract(frozenBackgroundVV);
var diffVH = Inew_linearVH.subtract(frozenBackgroundVH);

// Calculate absolute difference
var diffVVabs = diffVV.abs();
var diffVHabs = diffVH.abs();

// Calculate relative change (normalized by frozen background)
var relativeChangeVV = diffVVabs.divide(frozenBackgroundVV.max(0.01));
var relativeChangeVH = diffVHabs.divide(frozenBackgroundVH.max(0.01));

// Create change masks based on type
var changeVV, changeVH;

if (changeType === 'increase') {
  // Only detect increases (new > frozen)
  changeVV = diffVV.gt(0).and(relativeChangeVV.gt(changeThresholdVV));
  changeVH = diffVH.gt(0).and(relativeChangeVH.gt(changeThresholdVH));
} else if (changeType === 'decrease') {
  // Only detect decreases (new < frozen)
  changeVV = diffVV.lt(0).and(relativeChangeVV.gt(changeThresholdVV));
  changeVH = diffVH.lt(0).and(relativeChangeVH.gt(changeThresholdVH));
} else {
  // Detect both increases and decreases
  changeVV = relativeChangeVV.gt(changeThresholdVV);
  changeVH = relativeChangeVH.gt(changeThresholdVH);
}

// Combined change mask (change in either polarization)
var changeMask = changeVV.or(changeVH);

// Optional: Apply morphological operations to remove small patches
changeMask = changeMask.selfMask()
  .connectedPixelCount(100, true)
  .gte(minChangeArea);

// Calculate change magnitude for visualization
var changeMagnitude = relativeChangeVV.add(relativeChangeVH).divide(2);

// Calculate change direction (-1 to 1, negative = decrease, positive = increase)
var changeDirection = diffVV.add(diffVH).divide(2)
  .divide(frozenBackgroundVV.add(frozenBackgroundVH).divide(2).max(0.01));

// Print change statistics
var changeStats = changeMask.multiply(100).reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: geometry,
  scale: 10,
  maxPixels: 1e9
});
print('📊 Percentage of pixels with change:', changeStats);

// =========================================================================
// 5. Display Results
// =========================================================================

Map.setOptions('satellite');
var visParams = {min: [0, 0, 0], max: [1, 1, 1], gamma: 1};

// Display original layers
Map.addLayer(frozenBackgroundVV, {min: 0, max: 1}, 'Frozen Background VV', false);
Map.addLayer(meanVV, {min: 0, max: 1}, 'Mean VV', false);
Map.addLayer(frozenBackgroundVH, {min: 0, max: 1}, 'Frozen Background VH', false);
Map.addLayer(meanVH, {min: 0, max: 1}, 'Mean VH', false);

// Display change detection layers
Map.addLayer(relativeChangeVV, {min: 0, max: 0.5, palette: ['blue', 'yellow', 'red']}, 'Relative Change VV', false);
Map.addLayer(relativeChangeVH, {min: 0, max: 0.5, palette: ['blue', 'yellow', 'red']}, 'Relative Change VH', false);
Map.addLayer(changeMask.selfMask(), {palette: ['red']}, 'Change Mask', false);

// Display change direction (increase = red, decrease = blue)
var changeDirectionVis = changeDirection.updateMask(changeMask);
Map.addLayer(changeDirectionVis, {min: -0.5, max: 0.5, palette: ['blue', 'white', 'red']}, 'Change Direction', false);

// Create RGB composites
var rgbvv = ee.Image.cat(
  Inew_linearVV.max(frozenBackgroundVV.select([0])),
  frozenBackgroundVV.select([0]),
  frozenBackgroundVV.select([0])
).clip(simplifiedDonbas);

var rgbvh = ee.Image.cat(
  frozenBackgroundVH.select([0]).multiply(3),
  frozenBackgroundVH.select([0]).multiply(3),
  Inew_linearVH.max(frozenBackgroundVH.select([0])).multiply(3)
).clip(simplifiedDonbas);

// Original FBRpolar (full image)
var FBRpolar = ee.Image.cat(
  frozenBackgroundVV.select([0]),
  frozenBackgroundVH.select([0]).multiply(3),
  frozenBackgroundVH.select([0]).divide(frozenBackgroundVV.select([0]).max(0.0001))
).clip(simplifiedDonbas);

// FBRpolar with change mask applied (only change pixels visible)
var FBRpolarMasked = FBRpolar.updateMask(changeMask);

// Display both versions
Map.addLayer(FBRpolar, visParams, 'FBG in polarimetric colors (full)', false);
Map.addLayer(FBRpolarMasked, visParams, 'FBG in polarimetric colors (changes only)');

// Optional: Create a composite showing changes with context
// Changes in color, no-change areas in grayscale
var backgroundGray = frozenBackgroundVV.select([0]).multiply(0.5);
var contextComposite = FBRpolar.where(changeMask.not(), 
  ee.Image.cat(backgroundGray, backgroundGray, backgroundGray));
Map.addLayer(contextComposite, visParams, 'Changes with context', false);

// Create a special visualization for change types
// Red = increase, Blue = decrease, Intensity = magnitude
var increaseVis = diffVV.add(diffVH).divide(2).gt(0).and(changeMask);
var decreaseVis = diffVV.add(diffVH).divide(2).lt(0).and(changeMask);
var changeTypeVis = ee.Image.cat(
  increaseVis.multiply(changeMagnitude),  // Red channel
  ee.Image(0),                             // Green channel
  decreaseVis.multiply(changeMagnitude)    // Blue channel
);
Map.addLayer(changeTypeVis.updateMask(changeMask), {min: 0, max: 0.3}, 'Change Type Vis (R=increase, B=decrease)', false);

// =========================================================================
// 6. Change Detection Parameters UI
// =========================================================================

// Create a control panel for threshold adjustment
var thresholdPanel = ui.Panel({
  widgets: [
    ui.Label('🎛️ Change Detection Controls', {fontWeight: 'bold'}),
    ui.Label('Change Type:'),
    ui.Select({
      items: ['both', 'increase', 'decrease'],
      value: changeType,
      onChange: function(value) {
        changeType = value;
        updateChangeMask();
      }
    }),
    ui.Label('VV Threshold: ' + changeThresholdVV),
    ui.Slider({
      min: 0,
      max: 0.5,
      value: changeThresholdVV,
      step: 0.01,
      onChange: function(value) {
        changeThresholdVV = value;
        updateChangeMask();
      }
    }),
    ui.Label('VH Threshold: ' + changeThresholdVH),
    ui.Slider({
      min: 0,
      max: 0.5,
      value: changeThresholdVH,
      step: 0.01,
      onChange: function(value) {
        changeThresholdVH = value;
        updateChangeMask();
      }
    }),
    ui.Label('Min Area (pixels): ' + minChangeArea),
    ui.Slider({
      min: 1,
      max: 100,
      value: minChangeArea,
      step: 1,
      onChange: function(value) {
        minChangeArea = value;
        updateChangeMask();
      }
    })
  ],
  style: {
    position: 'top-left',
    padding: '8px'
  }
});

// Function to update change mask with new thresholds
function updateChangeMask() {
  // Recalculate differences
  var diffVV = Inew_linearVV.subtract(frozenBackgroundVV);
  var diffVH = Inew_linearVH.subtract(frozenBackgroundVH);
  var diffVVabs = diffVV.abs();
  var diffVHabs = diffVH.abs();
  var relativeChangeVV = diffVVabs.divide(frozenBackgroundVV.max(0.01));
  var relativeChangeVH = diffVHabs.divide(frozenBackgroundVH.max(0.01));
  
  // Apply change type filter
  if (changeType === 'increase') {
    changeVV = diffVV.gt(0).and(relativeChangeVV.gt(changeThresholdVV));
    changeVH = diffVH.gt(0).and(relativeChangeVH.gt(changeThresholdVH));
  } else if (changeType === 'decrease') {
    changeVV = diffVV.lt(0).and(relativeChangeVV.gt(changeThresholdVV));
    changeVH = diffVH.lt(0).and(relativeChangeVH.gt(changeThresholdVH));
  } else {
    changeVV = relativeChangeVV.gt(changeThresholdVV);
    changeVH = relativeChangeVH.gt(changeThresholdVH);
  }
  
  changeMask = changeVV.or(changeVH);
  
  // Apply size filter
  changeMask = changeMask.selfMask()
    .connectedPixelCount(100, true)
    .gte(minChangeArea);
  
  // Update masked layer
  FBRpolarMasked = FBRpolar.updateMask(changeMask);
  
  // Remove old layer and add new one
  Map.layers().forEach(function(layer) {
    if (layer.getName() === 'FBG in polarimetric colors (changes only)') {
      Map.remove(layer);
    }
  });
  
  Map.addLayer(FBRpolarMasked, visParams, 'FBG in polarimetric colors (changes only)');
}

Map.add(thresholdPanel);

// =========================================================================
// 7. Export Options
// =========================================================================

// Get the date string for file naming
var dateString = startDate.format('YYYY-MM-dd').getInfo() + '_to_' + endDate.format('YYYY-MM-dd').getInfo();

print('📤 Creating frozen background export tasks...');

// Export Frozen Background VV
Export.image.toDrive({
  image: frozenBackgroundVV.multiply(10000).toInt16(),
  description: 'FrozenBackground_VV_' + dateString,
  folder: 'GEE_exports',
  fileNamePrefix: 'FrozenBackground_VV_' + dateString,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export Frozen Background VH
Export.image.toDrive({
  image: frozenBackgroundVH.multiply(10000).toInt16(),
  description: 'FrozenBackground_VH_' + dateString,
  folder: 'GEE_exports',
  fileNamePrefix: 'FrozenBackground_VH_' + dateString,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export Change Mask
Export.image.toDrive({
  image: changeMask.toByte(),
  description: 'ChangeMask_FrozenBG_' + dateString,
  folder: 'GEE_exports',
  fileNamePrefix: 'ChangeMask_FrozenBG_' + dateString,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export RGB composite (full scene)
Export.image.toDrive({
  image: FBRpolar.multiply(255).toByte(),
  description: 'FBG_RGB_Full_' + dateString,
  folder: 'GEE_exports',
  fileNamePrefix: 'FBG_RGB_Full_' + dateString,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export RGB composite (changes only)
Export.image.toDrive({
  image: FBRpolarMasked.multiply(255).toByte(),
  description: 'FBG_RGB_ChangesOnly_' + dateString,
  folder: 'GEE_exports',
  fileNamePrefix: 'FBG_RGB_ChangesOnly_' + dateString,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export Change Direction
Export.image.toDrive({
  image: changeDirection.multiply(changeMask).multiply(10000).toInt16(),
  description: 'ChangeDirection_FrozenBG_' + dateString,
  folder: 'GEE_exports',
  fileNamePrefix: 'ChangeDirection_FrozenBG_' + dateString,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

print('📤 Export tasks created. Check the Tasks tab to run exports.');

// =========================================================================
// 8. Usage Instructions
// =========================================================================
print('');
print('🎛️ FROZEN BACKGROUND DETECTION PARAMETERS:');
print('• mu_threshold: CV threshold for frozen background computation (default: 0.2286)');
print('• durationMonths: Duration of time series in months (default: 12)');
print('• changeThresholdVV/VH: Thresholds for change detection (0-0.5)');
print('• minChangeArea: Minimum connected pixels to consider as change');
print('• changeType: Type of changes to detect (increase/decrease/both)');
print('');
print('💡 FROZEN BACKGROUND TIPS:');
print('• Frozen background removes temporal anomalies using CV threshold');
print('• Lower mu_threshold = more conservative background (removes more values)');
print('• Higher mu_threshold = more inclusive background (keeps more values)');
print('• Use UI controls to interactively adjust change detection parameters');
print('• Red areas = increases, Blue areas = decreases in change direction layer');
print('• Experiment with changeType to focus on specific change directions');
print('• Check frozen background layers to understand baseline conditions');
print('');
print('🔬 METHODOLOGY:');
print('• Based on Taillade et al. (2020) frozen background approach');
print('• Computes temporally stable background by removing CV outliers');
print('• Detects ephemeral changes by comparing new images to frozen background');
print('• Processes VV and VH polarizations separately for better sensitivity');
print('• Uses simplified statistical approach to avoid complex array operations');