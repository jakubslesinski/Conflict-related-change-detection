// =========================================================================
// 1. User parameters 
// =========================================================================

var str1='2022-02-24';     // Start of the Observation 'YY-MM-dd'
var str2='2025-06-30';     // End of the Observation 'YY-MM-dd'
var croppalet=0.66;

// THRESHOLDING PARAMETERS
var enableThresholding = true;        // Enable/disable thresholding
var cvThreshold = 0.25;              // Threshold for coefficient of variation (0.1-0.5)
var intensityThreshold = 0.02;       // Threshold for maximum intensity
var eigenThreshold = 0.3;            // Threshold for maximum eigenvalue (0.1-0.5)
var combinedThreshold = true;        // Use combination of all thresholds (true/false)
var useEigenThreshold = true;        // Whether to use eigenvalue threshold (true/false)

// -------------------------------------------------------------------------
// LOAD DONBAS REGION (you need to import your shapefile as an asset first)
// -------------------------------------------------------------------------

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

// Study area: Donbas region
var geometry = simplifiedDonbas.geometry();

// Center the map on the Donbas region
Map.centerObject(simplifiedDonbas, 13);
//Map.addLayer(simplifiedDonbas, {color: 'red'}, 'Donbas Region', false);
Map.addLayer(simplifiedDonbas.style({color: 'red', fillColor: '00000000'}));

// =========================================================================
// 2. COMPUTE COLLECTION
// =========================================================================

var startDate = ee.Date(str1);
var endDate = ee.Date(str2);
var ds = endDate.difference(startDate, 'day');

// Load the Sentinel-1 ImageCollection filtered by Donbas region
var sentinel1_liste = ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT')
  .filterDate(startDate, endDate)
  .filterBounds(simplifiedDonbas);

// Get orbit information
var NbOrbit = sentinel1_liste.aggregate_count_distinct('relativeOrbitNumber_start');
var ListOrbits = sentinel1_liste.aggregate_array('relativeOrbitNumber_start');
var freq = ee.Dictionary(ee.List(ListOrbits).reduce(ee.Reducer.frequencyHistogram()));
var array = ee.Array([freq.keys().map(ee.Number.parse), freq.values()]);
var frequences = array.slice(0,-1);
var arraysort = array.sort(frequences);
var index = ee.Number(NbOrbit).add(-1);
var orbite = arraysort.get([0,ee.Number(index)]);
var sentinel1 = sentinel1_liste.filterMetadata('relativeOrbitNumber_start', 'equals', orbite);

// Get metadata
var modePass = sentinel1.aggregate_array('orbitProperties_pass').distinct();
var acquisitionMode = sentinel1.aggregate_array('instrumentMode').distinct();

var allPolarizations = sentinel1.aggregate_array('transmitterReceiverPolarisation');
var commonPolarizations = allPolarizations.iterate(function(list, intersection) {
  return ee.List(intersection).removeAll(ee.List(intersection).removeAll(ee.List(list)));
}, allPolarizations.get(0));
var polarList = ee.List(commonPolarizations);

var NumberOfChannels = polarList.length();

// Check polarization requirements
if (NumberOfChannels.getInfo() !== 2) {
  throw new Error('⚠️ Polarimetric version not applicable: Exactly two polarimetric channels are required.');
}

var polar1 = polarList.get(0);
var polar2 = polarList.get(1);
var polar1 = ee.List([polar1]); 
var polar2 = ee.List([polar2]); 

// Checking and converting lists into strings
var formattedModePass = ee.String(ee.List(modePass).join(", "));
var formattedAcquisitionMode = ee.String(ee.List(acquisitionMode).join(", "));

// Check if `commonPolarizations` is empty and format the display
var formattedPolarizations = ee.Algorithms.If(
  ee.List(commonPolarizations).size().gt(0),  // Check if the list has elements
  ee.String(ee.List(commonPolarizations).join(", ")), 
  ee.String("None")  // If the list is empty, display "None".
);

var summary = ee.Dictionary({
  '🌍 Region': 'Donbas',
  '🛰️ Selected Orbit': orbite,
  '📡 Orbit Mode': formattedModePass,
  '🛠️ Acquisition Mode': formattedAcquisitionMode,
  '📊 Number of Polarimetric Channels': NumberOfChannels,
  '🎭 Common Polarizations in All Images': formattedPolarizations,
  '📅 Date Range': str1 + ' to ' + str2,
  '📸 Number of Images': sentinel1.size(),
  '🎚️ CV Threshold': cvThreshold,
  '📊 Intensity Threshold': intensityThreshold,
  '🔶 Eigenvalue Threshold': eigenThreshold,
  '🔧 Thresholding Enabled': enableThresholding,
  '👻 Transparency Mode': 'ON - non-change areas transparent'
});

print('📊 Sentinel-1 Polarimetric Summary:', summary);

// =========================================================================
// 3. Compute Multivariate Coefficient of Variation Min and Max bounds
// =========================================================================

// This function applies a linear scale conversion to each image in the collection
var amplitude = function(image) {
  var imlin = image.expression(
    'sqrt(intensity)', {
      'intensity': image
  });
  return imlin; // Conversion to linear scale for classical mean computation
};

// Selection of time series for polarimetric channels polar1 and polar2
var SeriesPolar1 = sentinel1.select(polar1).map(amplitude);
var SeriesPolar2 = sentinel1.select(polar2).map(amplitude);

// Computation of the mean for each polarimetric channel of time series 
var meanPolar1 = SeriesPolar1.mean();
var meanPolar2 = SeriesPolar2.mean();

// Normalization of time series by their respective means
var seriesPolar1 = SeriesPolar1.map(function(image) {
  return image.divide(meanPolar1);
});
var seriesPolar2 = SeriesPolar2.map(function(image) {
  return image.divide(meanPolar2);
});
var meanPolar1 = seriesPolar1.mean();
var meanPolar2 = seriesPolar2.mean();

var stdLinear1 = seriesPolar1.reduce(ee.Reducer.stdDev());
var magic1 = stdLinear1.divide(meanPolar1);
var stdLinear2 = seriesPolar2.reduce(ee.Reducer.stdDev());
var magic2 = stdLinear2.divide(meanPolar2);

// Computation of deviations from the mean for each time series
var deviationPolar1 = seriesPolar1.map(function(image) {
  return image.subtract(meanPolar1);
});

var deviationPolar2 = seriesPolar2.map(function(image) {
  return image.subtract(meanPolar2);
});

// Computation of c11 = mean((seriesPolar1 - meanPolar1)^2)
var c11 = deviationPolar1.map(function(image) {
  return image.pow(2);  // (seriesPolar1 - meanPolar1)^2
}).reduce(ee.Reducer.mean());

// Manual computation of c22 = mean((seriesPolar2 - meanPolar2)^2)
var c22 = deviationPolar2.map(function(image) {
  return image.pow(2);  // (seriesPolar2 - meanPolar2)^2
}).reduce(ee.Reducer.mean());

// Pair the two image series using zip()
var pairedSeries = deviationPolar1.toList(deviationPolar1.size()).zip(deviationPolar2.toList(deviationPolar2.size()));

// Function to calculate pixel-wise product for each pair of images
var productSeries = ee.ImageCollection(pairedSeries.map(function(pair) {
  pair = ee.List(pair);
  var image1 = ee.Image(pair.get(0));
  var image2 = ee.Image(pair.get(1));
  return image1.multiply(image2).copyProperties(image1, image1.propertyNames());
}));
var c12 = productSeries.mean();

// Computation of delta = sqrt((c11 - c22)^2 + 4 * c12^2)
var delta = c11.subtract(c22).pow(2).add(c12.pow(2).multiply(4)).sqrt();

// Computation of eigenvalues (maximum and minimum)
var lambdamaxi = c11.add(c22).add(delta).divide(2);  // Maximum eigenvalue
var lambdamini = c11.add(c22).subtract(delta).divide(2);  // Minimum eigenvalue

// Computation of the norm of the mean vector
var normMU = (meanPolar1.pow(2).add(meanPolar2.pow(2))).sqrt();

// Normalization of eigenvalues by the norm of the mean vector
var limitmin_per_pixel = lambdamini.sqrt().divide(normMU);
var limitmax_per_pixel = lambdamaxi.sqrt().divide(normMU);

var magic = magic1.max(magic2);

// Clip all results to Donbas region
limitmin_per_pixel = limitmin_per_pixel.clip(simplifiedDonbas);
limitmax_per_pixel = limitmax_per_pixel.clip(simplifiedDonbas);
magic1 = magic1.clip(simplifiedDonbas);
magic2 = magic2.clip(simplifiedDonbas);
magic = magic.clip(simplifiedDonbas);

// =========================================================================
// 4. Compute CLASSICAL REACTIV for Comparison
// =========================================================================

var imax1 = sentinel1.select(polar1).max();
var imax2 = sentinel1.select(polar2).max();
var imax = imax1.max(imax2);

// This function affects value of days for pixels where maximum is reached
var time1 = function(image) {
  var days = image.date().difference(startDate, 'day').divide(ds); //divide by the period of time observed
  return image.where(image.lt(imax1),0).where(image.gte(imax1),days);
};
var days1 = sentinel1.select(polar1).map(time1).sum();

var time2 = function(image) {
  var days = image.date().difference(startDate, 'day').divide(ds); //divide by the period of time observed
  return image.where(image.lt(imax2),0).where(image.gte(imax2),days);
};
var days2 = sentinel1.select(polar2).map(time2).sum();
var days = (days2.where(magic2.lte(magic1),days1));

// Images of Number of images: sizepile
var unit = function(image) {
  var imunit = image.multiply(0).add(1);
  return imunit; // conversion in linear, then compute mean: classical mean
};
var sizepile = sentinel1.select(polar1).map(unit).sum(); 

// Parameter for dynamics
var mu = 0.2286; // Theoretical mean for Rayleigh Nakagam L=4.9
var stdmu = ee.Image(0.1616);
stdmu = stdmu.divide(sizepile.sqrt()); // Theoretical std for Rayleigh Nakagami L=4.9
var magicnorm = magic.subtract(mu).divide(stdmu.multiply(10)).clamp(0,1);

var rgb = ee.Image.cat(days.multiply(croppalet),magicnorm,imax.clamp(0,1)).hsvToRgb();

// Clip final results
imax = imax.clip(simplifiedDonbas);
days = days.clip(simplifiedDonbas);
magicnorm = magicnorm.clip(simplifiedDonbas);
rgb = rgb.clip(simplifiedDonbas);

// =========================================================================
// 5. POLARIMETRIC THRESHOLDING - NEW FUNCTIONALITY
// =========================================================================

// Creating thresholding masks
var cvMask = magic.gte(cvThreshold);
var intensityMask = imax.gte(intensityThreshold);
var eigenMask = limitmax_per_pixel.gte(eigenThreshold);

// Combining masks depending on settings
var changeMask;
if (combinedThreshold) {
  if (useEigenThreshold) {
    changeMask = cvMask.and(intensityMask).and(eigenMask);
    print('🎭 Using combined CV, intensity, and eigenvalue thresholding');
  } else {
    changeMask = cvMask.and(intensityMask);
    print('🎭 Using combined CV and intensity thresholding');
  }
} else {
  changeMask = cvMask;
  print('🎭 Using CV thresholding only');
}

// Creating visualization with thresholding
var rgbThresholded;
var rgbImage3Thresholded;

if (enableThresholding) {
  // Apply transparency mask - only changes will be visible
  rgbThresholded = rgb.updateMask(changeMask);
  
  // Apply thresholding to lmin lmax CVmax composite
  var rgbImage3 = ee.Image.cat(limitmin_per_pixel.clamp(0,0.25).divide(0.25), 
                              limitmax_per_pixel.clamp(0,0.5).divide(0.5), 
                              magic.clamp(0,0.75).divide(0.75)).clip(simplifiedDonbas);
  rgbImage3Thresholded = rgbImage3.updateMask(changeMask);
  
  // Masking statistics
  var totalPixels = changeMask.select([0]).reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: simplifiedDonbas.geometry(),
    scale: 100,
    maxPixels: 1e9
  });
  
  var changePixels = changeMask.select([0]).reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: simplifiedDonbas.geometry(),
    scale: 100,
    maxPixels: 1e9
  });
  
  print('📈 Polarimetric Thresholding Statistics:');
  print('Total pixels:', totalPixels);
  print('Change pixels:', changePixels);
  print('🔍 Areas without changes will be transparent');
  
} else {
  rgbThresholded = rgb;
  var rgbImage3 = ee.Image.cat(limitmin_per_pixel.clamp(0,0.25).divide(0.25), 
                              limitmax_per_pixel.clamp(0,0.5).divide(0.5), 
                              magic.clamp(0,0.75).divide(0.75)).clip(simplifiedDonbas);
  rgbImage3Thresholded = rgbImage3;
  print('🚫 Thresholding disabled - showing all data');
}

// =========================================================================
// 6. VISUALIZATION
// =========================================================================

Map.setOptions('satellite');
var visparams = {min: [0, 0, 0],max: [1, 1, 1],gamma: 2};

var couchemin = limitmin_per_pixel.clamp(0,0.25).divide(0.25).pow(1);
var couchemax = limitmax_per_pixel.clamp(0,0.5).divide(0.5).pow(1);

// Layers without thresholding (disabled by default)
Map.addLayer(couchemin, {min: 0, max: 1, palette: ['black', 'white']}, 'Limit Min', false);
Map.addLayer(couchemax, {min: 0, max: 1, palette: ['black', 'white']}, 'Limit Max', false);
Map.addLayer(magic1.clamp(0,0.5).divide(0.5).pow(2), {min: 0, max: 1, gamma: 1}, 'CV, Polar 1', false);
Map.addLayer(magic2.clamp(0,0.5).divide(0.5).pow(2), {min: 0, max: 1, gamma: 1}, 'CV, Polar 2', false);

// Main layers
Map.addLayer(rgb, visparams,'REACTIV Original', false);
Map.addLayer(rgbThresholded, visparams,'REACTIV Transparent');

var rgbImage3 = ee.Image.cat(couchemin, couchemax, magic.clamp(0,0.75).divide(0.75)).clip(simplifiedDonbas);
Map.addLayer(rgbImage3, {min: [0, 0, 0],max: [1, 1, 1],gamma: 0.5}, 'lmin lmax CVmax Original', false);
Map.addLayer(rgbImage3Thresholded, {min: [0, 0, 0],max: [1, 1, 1],gamma: 0.5}, 'lmin lmax CVmax Transparent');

// Add masks as additional layers for control
Map.addLayer(changeMask.selfMask(), {palette: ['white']}, 'Change Mask', false);
Map.addLayer(cvMask.selfMask(), {palette: ['yellow']}, 'CV Mask', false);
Map.addLayer(intensityMask.selfMask(), {palette: ['cyan']}, 'Intensity Mask', false);
if (useEigenThreshold) {
  Map.addLayer(eigenMask.selfMask(), {palette: ['magenta']}, 'Eigenvalue Mask', false);
}

// =========================================================================
// 7. Export Options
// =========================================================================

// Get the date string for file naming
var dateString = startDate.format('YYYY-MM-dd').getInfo() + '_to_' + endDate.format('YYYY-MM-dd').getInfo();
var thresholdSuffix = enableThresholding ? '_thresh_' + cvThreshold.toString().replace('.', '') : '';

// Export REACTIV RGB visualization (original)
Export.image.toDrive({
  image: rgb.multiply(255).toByte(),
  description: 'REACTIV_RGB_' + dateString,
  folder: 'GEE_exports',
  fileNamePrefix: 'REACTIV_RGB_' + dateString,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export REACTIV RGB visualization (with thresholding and transparency)
Export.image.toDrive({
  image: rgbThresholded.multiply(255).toByte(),
  description: 'REACTIV_RGB_Transparent_' + dateString + thresholdSuffix,
  folder: 'GEE_exports',
  fileNamePrefix: 'REACTIV_RGB_Transparent_' + dateString + thresholdSuffix,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export lmin lmax CVmax RGB (original)
Export.image.toDrive({
  image: rgbImage3.multiply(255).toByte(),
  description: 'LminLmaxCVmax_RGB_' + dateString,
  folder: 'GEE_exports',
  fileNamePrefix: 'LminLmaxCVmax_RGB_' + dateString,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export lmin lmax CVmax RGB (with thresholding and transparency)
Export.image.toDrive({
  image: rgbImage3Thresholded.multiply(255).toByte(),
  description: 'LminLmaxCVmax_RGB_Transparent_' + dateString + thresholdSuffix,
  folder: 'GEE_exports',
  fileNamePrefix: 'LminLmaxCVmax_RGB_Transparent_' + dateString + thresholdSuffix,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export change mask
Export.image.toDrive({
  image: changeMask.toByte(),
  description: 'Change_Mask_Polarimetric_' + dateString + thresholdSuffix,
  folder: 'GEE_exports',
  fileNamePrefix: 'Change_Mask_Polarimetric_' + dateString + thresholdSuffix,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export Limit Min
Export.image.toDrive({
  image: limitmin_per_pixel.multiply(10000).toInt16(),
  description: 'LimitMin_' + dateString,
  folder: 'GEE_exports',
  fileNamePrefix: 'LimitMin_' + dateString,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export Limit Max
Export.image.toDrive({
  image: limitmax_per_pixel.multiply(10000).toInt16(),
  description: 'LimitMax_' + dateString,
  folder: 'GEE_exports',
  fileNamePrefix: 'LimitMax_' + dateString,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export CV Polar 1
Export.image.toDrive({
  image: magic1.multiply(10000).toInt16(),
  description: 'CV_Polar1_' + dateString,
  folder: 'GEE_exports',
  fileNamePrefix: 'CV_Polar1_' + dateString,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export CV Polar 2
Export.image.toDrive({
  image: magic2.multiply(10000).toInt16(),
  description: 'CV_Polar2_' + dateString,
  folder: 'GEE_exports',
  fileNamePrefix: 'CV_Polar2_' + dateString,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export days (temporal information)
Export.image.toDrive({
  image: days.multiply(10000).toInt16(),
  description: 'Days_MaxIntensity_' + dateString,
  folder: 'GEE_exports',
  fileNamePrefix: 'Days_MaxIntensity_' + dateString,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export maximum intensity
Export.image.toDrive({
  image: imax.multiply(10000).toInt16(),
  description: 'MaxIntensity_' + dateString,
  folder: 'GEE_exports',
  fileNamePrefix: 'MaxIntensity_' + dateString,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

print('📤 Export tasks created. Check the Tasks tab to run exports.');

// =========================================================================
// 8. TIME PROFILE PLOT
// =========================================================================

// Create the title label.
var title = ui.Label('Click to inspect (Polarimetric Thresholding: ' + (enableThresholding ? 'ON' : 'OFF') + ')');
title.style().set('position', 'bottom-center');
Map.add(title);

// Create a panel to hold the chart.
var panel = ui.Panel();
panel.style().set({
  width: '400px',
  position: 'bottom-right'
});
Map.add(panel);

var selectedBands = ee.List([
  ee.Algorithms.If(polar1, polar1.get(0), null), 
  ee.Algorithms.If(polar2, polar2.get(0), null)
]).removeAll([null]);

// Convert `selectedBands` to a native JavaScript array
var selectedBandsJS = selectedBands.getInfo();  

// Register a function to draw a chart when a user clicks on the map.
Map.onClick(function(coords) {
  panel.clear();
  var point = ee.Geometry.Point(coords.lon, coords.lat);
  
  // Check if the clicked point is within Donbas region
  var withinDonbas = simplifiedDonbas.geometry().contains(point);
  var isInside = withinDonbas.getInfo();
  
  if (!isInside) {
    panel.add(ui.Label('⚠️ Click within the Donbas region to see temporal profiles'));
    return;
  }
  
  // Get values at point
  var pointValues = ee.Image.cat(magic, imax, limitmax_per_pixel, changeMask).reduceRegion({
    reducer: ee.Reducer.first(),
    geometry: point,
    scale: 10
  }).getInfo();
  
  var cvValue = pointValues.cv || 0;
  var intensityValue = pointValues.intensity || 0;
  var eigenValue = pointValues.limitmax_per_pixel || 0;
  var isChange = pointValues.cv_1 || 0;
  
  var infoLabel = ui.Label('CV: ' + cvValue.toFixed(3) + 
                          ' | Intensity: ' + intensityValue.toFixed(3) + 
                          ' | Eigenvalue: ' + eigenValue.toFixed(3) +
                          ' | Change detected: ' + (isChange ? 'YES' : 'NO'));
  panel.add(infoLabel);
  
  print('Coordinates for the plot: ', point);
  
  if (selectedBandsJS.length === 0) {
    if (panel.widgets().length() === 0) { // Prevent multiple display of message
      panel.add(ui.Label("❌ Error: No valid band for graph display."));
    }
    return; // Stop Execution Here
  }
  
  // Dynamic graph creation based on available bands
  var chart2 = ui.Chart.image.series(sentinel1.select(selectedBandsJS), point, null, 30)
    .setChartType('ScatterChart')
    .setOptions({
      title: 'Temporal Profile at: ' + coords.lon.toFixed(4) + ', ' + coords.lat.toFixed(4),
      hAxis: { title: 'Acquisition Date' },
      vAxis: { title: 'Intensity Values (linear)' },
    });

  panel.add(chart2);
});
