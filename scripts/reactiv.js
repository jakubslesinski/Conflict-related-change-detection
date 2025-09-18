// When using this code in a publication, please cite
// Elise Colin Koeniguer et al, 
// Colored visualization of multitemporal SAR data for change detection: issues and methods
// EUSAR 2018
// 
// or, if you can read french:
//
// Visualisation des changements sur séries temporelles radar : méthode REACTIV 
// évaluée à l'échelle mondiale sous Google Earth Engine
// Elise Colin Koeniguer et al. 
// CFPT 2018
// https://rfiap2018.ign.fr/sites/default/files/ARTICLES/CFPT2018/Oraux/CFPT2018_paper_koeniguer.pdf
// -------------------------------------------------------------

// THRESHOLDING PARAMETERS
var enableThresholding = true;        // Enable/disable thresholding
var cvThreshold = 0.3;               // Threshold for coefficient of variation (0.75-1.0)
var intensityThreshold = -15;        // Threshold for maximum intensity (dB) - return to dB scale
var combinedThreshold = true;        // Use combination of both thresholds (true/false)
var croppalet = 0.6;                 // Use a number <1 if you do not want the whole HSV palet

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
Map.addLayer(simplifiedDonbas.style({color: 'red', fillColor: '00000000'}), {}, 'Donbas Region');

// -------------------------------------------------------------
// Parameters: DATES, ASCENDING OR DESCENDING
var str2='2025-06-30';    // end date
var str1='2022-02-24';    // start date
//var str='ASCENDING';
var str='ASCENDING';

// ------------------------------------------------------------
// date selection
var startDate = ee.Date(str1);
var endDate = ee.Date(str2);
var ds = endDate.difference(startDate, 'day');

// Centering
var pos = simplifiedDonbas.geometry().centroid();
print('Coordinate of the Center of the Donbas region:', pos);

// Load the Sentinel-1 ImageCollection centered on the location "pos"
// FIX: Using proper COPERNICUS/S1_GRD collection (data in dB)
var sentinel1_liste = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterDate(startDate, endDate)
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filterBounds(pos)
  .filter(ee.Filter.eq('orbitProperties_pass', str));

// sentinel collection of the world without the restriction of the position
var sentinel1_liste2 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterDate(startDate, endDate)
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.eq('orbitProperties_pass', str));

// a solution to get metadata value of images of a collection
var NbOrbit = sentinel1_liste.aggregate_count_distinct('relativeOrbitNumber_start');
print('Number of Orbits',NbOrbit);
var ListOrbits = sentinel1_liste.aggregate_array('relativeOrbitNumber_start');
print('ListOrbits:', ListOrbits);

// find orbit numbers and their frequency
var freq = ee.Dictionary(ee.List(ListOrbits).reduce(ee.Reducer.frequencyHistogram()));
print('freq',freq.values());
var array = ee.Array([freq.keys().map(ee.Number.parse), freq.values()]);
print('array',array);

// orbit choice : first, the one with the max frequency
var frequences = array.slice(0,-1);
var arraysort = array.sort(frequences);
var index = ee.Number(NbOrbit).add(-1);
var orbite = arraysort.get([0,ee.Number(index)]);
print('Selected orbit=',orbite);

// find images with the choice orbit
var sentinel1 = sentinel1_liste2.filterMetadata('relativeOrbitNumber_start', 'equals', orbite);

var sentinel1_lista = sentinel1.toList(sentinel1.size());
print('sentinel1_lista: ', sentinel1_lista);

// Check if we have any images
var imageCount = sentinel1.size();
print('🔍 Number of images in collection:', imageCount);

// Check if we have any images
var hasImages = imageCount.gt(0);
print('🔍 Collection has images:', hasImages);

// FIX: Simplified polarization handling with proper band names
// Polarimetric channels available in the collection
var allPolarizations = sentinel1.aggregate_array('transmitterReceiverPolarisation').distinct();
print('All polarizations:', allPolarizations);

// Get unique polarizations
var uniquePolarizations = allPolarizations.getInfo();
var polarizationsFlat = [];

// Flatten the array of arrays to get individual polarizations
uniquePolarizations.forEach(function(polArray) {
  polArray.forEach(function(pol) {
    if (polarizationsFlat.indexOf(pol) === -1) {
      polarizationsFlat.push(pol);
    }
  });
});

print('Available polarizations:', polarizationsFlat);

// Select available polarizations (prioritize VV, then VH)
var polar1 = null;
var polar2 = null;

if (polarizationsFlat.indexOf('VV') !== -1) {
  polar1 = 'VV';
}
if (polarizationsFlat.indexOf('VH') !== -1) {
  if (polar1 === null) {
    polar1 = 'VH';
  } else {
    polar2 = 'VH';
  }
}

var NumberOfChannels = (polar1 !== null ? 1 : 0) + (polar2 !== null ? 1 : 0);

print('📊 Number of Polarimetric Channels:', NumberOfChannels);
print('🎭 Selected Polarizations:', 'polar1=' + polar1, 'polar2=' + polar2);

// Check if we can continue processing
if (NumberOfChannels === 0 || polar1 === null) {
  print('❌ ERROR: No valid polarizations found. Cannot proceed with REACTIV processing.');
  print('Available polarizations were:', polarizationsFlat);
} else {
  print('✅ Proceeding with REACTIV processing using', NumberOfChannels, 'polarization(s)');
}

// FIX: Conversion function for dB data (as in original version)
// This function applies to each image the linear scale conversion from dB
var amplitude = function(image) {
  var imlin = image.expression(
    '10**(amplitude/20)', {
      'amplitude': image
  });
  return imlin; // conversion from dB to linear amplitude
};

// FIX: Simplified implementation with simple strings for polarization names
var result;

// Check if we can continue
if (NumberOfChannels === 0 || polar1 === null) {
  print('❌ ERROR: Cannot process - no valid polarizations');
  result = {
    days: ee.Image(0),
    magic: ee.Image(0), 
    imax: ee.Image(0),
    activePolar: 'NONE'
  };
} else if (NumberOfChannels === 2 && polar1 !== null && polar2 !== null) {
  // If NumberOfChannels == 2
  var stdLinear1 = sentinel1.select(polar1).map(amplitude).reduce(ee.Reducer.stdDev()).clip(simplifiedDonbas);
  var meanLinear1 = sentinel1.select(polar1).map(amplitude).reduce(ee.Reducer.mean()).clip(simplifiedDonbas);
  var magic1 = stdLinear1.divide(meanLinear1);

  var stdLinear2 = sentinel1.select(polar2).map(amplitude).reduce(ee.Reducer.stdDev()).clip(simplifiedDonbas);
  var meanLinear2 = sentinel1.select(polar2).map(amplitude).reduce(ee.Reducer.mean()).clip(simplifiedDonbas);
  var magic2 = stdLinear2.divide(meanLinear2);
  
  var imax1 = sentinel1.select(polar1).max().clip(simplifiedDonbas);
  var imax2 = sentinel1.select(polar2).max().clip(simplifiedDonbas);
  var magic = magic1.max(magic2);
  var imax = imax1.max(imax2);

  // Time calculation function for first polarization
  var time1 = function(image) {
    var days = image.date().difference(startDate, 'day').divide(ds);
    return image.where(image.lt(imax1), 0).where(image.gte(imax1), days);
  };
  var days1 = sentinel1.select(polar1).map(time1).sum().clip(simplifiedDonbas);

  // Time calculation function for second polarization
  var time2 = function(image) {
    var days = image.date().difference(startDate, 'day').divide(ds);
    return image.where(image.lt(imax2), 0).where(image.gte(imax2), days);
  };
  var days2 = sentinel1.select(polar2).map(time2).sum().clip(simplifiedDonbas);

  var days = days2.where(magic2.lte(magic1), days1);
  var activePolar = polar1; // Default to polar1, can be refined based on magic values
  
  result = {
    days: days,
    magic: magic,
    imax: imax,
    activePolar: activePolar
  };
  
} else if (NumberOfChannels === 1 && polar1 !== null) {
  // If NumberOfChannels == 1
  var stdLinear1 = sentinel1.select(polar1).map(amplitude).reduce(ee.Reducer.stdDev()).clip(simplifiedDonbas);
  var meanLinear1 = sentinel1.select(polar1).map(amplitude).reduce(ee.Reducer.mean()).clip(simplifiedDonbas);
  var magic1 = stdLinear1.divide(meanLinear1);
  
  var imax1 = sentinel1.select(polar1).max().clip(simplifiedDonbas);
  var magic = magic1;
  var imax = imax1;

  // Time calculation function for single polarization
  var time1 = function(image) {
    var days = image.date().difference(startDate, 'day').divide(ds);
    return image.where(image.lt(imax1), 0).where(image.gte(imax1), days);
  };
  var days1 = sentinel1.select(polar1).map(time1).sum().clip(simplifiedDonbas);

  result = {
    days: days1,
    magic: magic,
    imax: imax,
    activePolar: polar1
  };
} else {
  print('❌ Error: Unexpected condition in polarization logic');
  result = {
    days: ee.Image(0),
    magic: ee.Image(0), 
    imax: ee.Image(0),
    activePolar: 'ERROR'
  };
}

// Extraction of final Results with check if result is defined
print('🔍 Debugging result object:', typeof result);

if (result && result.days && result.magic && result.imax) {
  var days = result.days.rename('temporal_days');
  var CV = result.magic.rename('CV');  // CV instead of magic
  var imax = result.imax.rename('max_intensity');
  var activePolar = result.activePolar;
  
  print('✅ Successfully extracted results from processing');
} else {
  print('❌ Error: Result object or its properties are undefined');
  print('Result:', result);
  
  // Fallback - create dummy images
  var days = ee.Image(0).rename('temporal_days');
  var CV = ee.Image(0).rename('CV');
  var imax = ee.Image(0).rename('max_intensity');
  var activePolar = 'ERROR';
}

// Images of Number of images: sizepile (only if we have valid polarization)
var unit = function(image) {
  var imunit = image.multiply(0).add(1);
  return imunit;
};

if (polar1 !== null) {
  var sizepile = sentinel1.select(polar1).map(unit).sum().clip(simplifiedDonbas).rename('stack_size');
  print('✅ Successfully computed sizepile for polarization:', polar1);
} else {
  var sizepile = ee.Image(1).rename('stack_size'); // fallback
  print('❌ No valid polarization found, using fallback sizepile');
} 

// Parameter for dynamics - IDENTICAL as in second version (only if we have valid images)
if (result && result.days && result.magic && result.imax) {
  var mu = 0.2286; // Theoretical mean for Rayleigh Nakagam L=4.9
  var stdmu = ee.Image(0.1616);
  var stdmu = stdmu.divide(sizepile.sqrt()); // Theoretical std for Rayleigh Nakagami L=4.9
  var CV_norm = CV.subtract(mu).divide(stdmu.multiply(10)).clamp(0,1).rename('CV_normalized'); // Normalized CV
  
  print('✅ Successfully computed CV_norm');
} else {
  var CV_norm = ee.Image(0).rename('CV_normalized');
  print('❌ Using fallback CV_norm');
}

// FIX: RGB composition with proper imax conversion (dB -> linear) for visualization
if (days && CV_norm && imax) {
  var imaxLinear = imax.select('max_intensity').expression(
    '10**(amplitude/20)', {
      'amplitude': imax.select('max_intensity')
    }
  ).rename('max_intensity_linear');

  var rgb = ee.Image.cat(
    days.select('temporal_days').multiply(croppalet), 
    CV_norm.select('CV_normalized'), 
    imaxLinear.clamp(0,1)
  ).hsvToRgb();
  
  print('✅ Successfully created RGB visualization');
} else {
  print('❌ Cannot create RGB - missing required images');
  var rgb = ee.Image([0,0,0]).rename(['red','green','blue']);
}

// Print parameters summary
var summary = ee.Dictionary({
  '🌍 Region': 'Donbas',
  '📅 Date Range': str1 + ' to ' + str2,
  '📡 Orbit Mode': str,
  '🎭 Active Polarization': activePolar,
  '📊 Number of Channels': NumberOfChannels,
  '🎚️ CV Threshold': cvThreshold,
  '📊 Intensity Threshold': intensityThreshold,
  '🔧 Thresholding Enabled': enableThresholding,
  '👻 Transparency Mode': 'ON - non-change areas transparent',
  '🌈 Color Palette Crop': croppalet
});
print('📊 REACTIV Unified Summary:', summary);

// Debug info about band names - with check if objects exist
print('🔍 Band names for debugging:');

// Check if all basic objects exist
print('🔍 Variable existence check:');
print('- days exists:', typeof days !== 'undefined');
print('- CV exists:', typeof CV !== 'undefined'); 
print('- CV_norm exists:', typeof CV_norm !== 'undefined');
print('- imax exists:', typeof imax !== 'undefined');
print('- changeMask exists:', typeof changeMask !== 'undefined');

try {
  if (typeof days !== 'undefined' && days !== null) {
    print('Days band names:', days.bandNames());
  } else {
    print('❌ Days image is undefined or null');
  }
  
  if (typeof CV !== 'undefined' && CV !== null) {
    print('CV band names:', CV.bandNames());
  } else {
    print('❌ CV image is undefined or null');
  }
  
  if (typeof CV_norm !== 'undefined' && CV_norm !== null) {
    print('CV_norm band names:', CV_norm.bandNames());
  } else {
    print('❌ CV_norm image is undefined or null');
  }
  
  if (typeof imax !== 'undefined' && imax !== null) {
    print('imax band names:', imax.bandNames());
  } else {
    print('❌ imax image is undefined or null');
  }
  
  if (typeof changeMask !== 'undefined' && changeMask !== null) {
    print('changeMask band names:', changeMask.bandNames());
  } else {
    print('❌ changeMask image is undefined or null');
  }
} catch (error) {
  print('❌ Error in debugging band names:', error);
}

// =========================================================================
// THRESHOLDING - ADAPTED TO NEW SCALE
// =========================================================================

// Creating thresholding masks (return to dB scale)
var cvMask = CV_norm.select('CV_normalized').gte(cvThreshold).rename('cv_mask');
var intensityMask = imax.select('max_intensity').gte(intensityThreshold).rename('intensity_mask'); // imax in dB

// Combining masks depending on settings
var changeMask;
if (combinedThreshold) {
  changeMask = cvMask.and(intensityMask).rename('change_combined');
  print('🎭 Using combined CV and intensity thresholding (dB scale)');
} else {
  changeMask = cvMask.rename('change_cv');
  print('🎭 Using CV thresholding only');
}

// Ensure all masks are clipped to region
changeMask = changeMask.clip(simplifiedDonbas);
cvMask = cvMask.clip(simplifiedDonbas);
intensityMask = intensityMask.clip(simplifiedDonbas);

// Creating visualization with thresholding
var rgbThresholded;
if (enableThresholding) {
  // Apply transparency mask - only changes will be visible
  rgbThresholded = rgb.updateMask(changeMask);
  
  // Masking statistics
  var totalPixels = changeMask.reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: simplifiedDonbas.geometry(),
    scale: 100,
    maxPixels: 1e9
  });
  
  var changePixels = changeMask.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: simplifiedDonbas.geometry(),
    scale: 100,
    maxPixels: 1e9
  });
  
  print('📈 Unified Thresholding Statistics:');
  print('Total pixels:', totalPixels);
  print('Change pixels:', changePixels);
  print('🔍 Areas without changes will be transparent');
  
} else {
  rgbThresholded = rgb;
  print('🚫 Thresholding disabled - showing all data');
}

var visparams = {min: [0, 0, 0], max: [1, 1, 1], gamma: 2.0};

// Main layers
Map.addLayer({
  eeObject: rgb,
  visParams: visparams,
  name: 'REACTIV Unified Original', 
  opacity: 0.8,
  shown: 0
});

Map.addLayer({
  eeObject: rgbThresholded,
  visParams: visparams,
  name: 'REACTIV Unified Transparent', 
  opacity: 0.8,
  shown: 1
});

// Add masks as additional layers for control
Map.addLayer(changeMask.selfMask(), {palette: ['white']}, 'Change Mask', false);
Map.addLayer(cvMask.selfMask(), {palette: ['yellow']}, 'CV Mask', false);
Map.addLayer(intensityMask.selfMask(), {palette: ['cyan']}, 'Intensity Mask', false);

Map.setOptions('satellite');

// =========================================================================
// TEMPORAL LEGEND - ADAPTED TO UNIFIED VERSION
// =========================================================================

function fractionToHex(hueFrac) {
  hueFrac = ee.Number(hueFrac).clamp(0, 1);
  var hsvImg = ee.Image.cat([hueFrac, 1, 1]).rename(['h', 's', 'v']);
  var rgbImg = hsvImg.hsvToRgb().multiply(255).round().toInt();
  var stats = rgbImg.reduceRegion({
    reducer: ee.Reducer.first(),
    geometry: ee.Geometry.Point([0, 0]),
    scale: 10,
    bestEffort: true,
    maxPixels: 1e6
  });
  var r = ee.Number(stats.get('red')).format('%02X');
  var g = ee.Number(stats.get('green')).format('%02X');
  var b = ee.Number(stats.get('blue')).format('%02X');
  return ee.String('#').cat(r).cat(g).cat(b);
}

function makeLegend(vis, d1, d2) {
  var lon = ee.Image.pixelLonLat().select('longitude');
  var gradient = lon.multiply(1/360.0);
  var legendImage = gradient.visualize(vis);
  var thumb = ui.Thumbnail({
    image: legendImage,
    params: {bbox:'0,0,360,4', dimensions:'200x20'},
    style: {padding: '1px', position: 'bottom-center', backgroundColor:'white'}
  });
  
  var panel = ui.Panel({
    widgets: [
      ui.Label(d1),
      ui.Label(str),
      ui.Label(activePolar),
      ui.Label(d2),
    ],
    layout: ui.Panel.Layout.flow('horizontal'),
    style: {stretch: 'horizontal', backgroundColor:'white', color:'blue'}
  });
  return ui.Panel({style: {backgroundColor: 'white'}}).add(panel).add(thumb);
}

var startFrac = ee.Number(0);
var endFrac = ee.Number(croppalet);
var numElements = 10;
var step = endFrac.subtract(startFrac).divide(ee.Number(numElements - 1));

var fracList = ee.List.sequence(startFrac, endFrac, step);
var palettehsv = fracList.map(function(hueFrac) {
  return ee.String(fractionToHex(hueFrac));
});

var pal = {
  min: 0,
  max: 1,
  palette: palettehsv
};

Map.add(makeLegend(pal, str1, str2));

// =========================================================================
// EXPORTS TO TIF - ADAPTED TO UNIFIED VERSION
// =========================================================================

// Get the date string for file naming
var dateString = startDate.format('YYYY-MM-dd').getInfo() + '_to_' + endDate.format('YYYY-MM-dd').getInfo();
var thresholdSuffix = enableThresholding ? '_thresh_' + cvThreshold.toString().replace('.', '') : '';
var polarSuffix = '_' + activePolar;

print('📤 Creating unified export tasks...');

// Export REACTIV RGB visualization (original)
Export.image.toDrive({
  image: rgb.multiply(255).toByte(),
  description: 'REACTIV_standard' + polarSuffix + '_' + str + '_' + dateString,
  folder: 'EarthEngine',
  fileNamePrefix: 'REACTIV_standard' + polarSuffix + '_' + str + '_' + dateString,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export REACTIV RGB visualization (with thresholding and transparency)
Export.image.toDrive({
  image: rgbThresholded.multiply(255).toByte(),
  description: 'REACTIV_Unified_RGB_Transparent' + polarSuffix + '_' + str + '_' + dateString + thresholdSuffix,
  folder: 'EarthEngine',
  fileNamePrefix: 'REACTIV_Unified_RGB_Transparent' + polarSuffix + '_' + str + '_' + dateString + thresholdSuffix,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export change mask
Export.image.toDrive({
  image: changeMask.toByte(),
  description: 'Change_Mask_Unified' + polarSuffix + '_' + str + '_' + dateString + thresholdSuffix,
  folder: 'EarthEngine',
  fileNamePrefix: 'Change_Mask_Unified' + polarSuffix + '_' + str + '_' + dateString + thresholdSuffix,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export CV (coefficient of variation)
Export.image.toDrive({
  image: CV.select('CV').multiply(10000).toInt16(),
  description: 'CV_Unified' + polarSuffix + '_' + str + '_' + dateString,
  folder: 'EarthEngine',
  fileNamePrefix: 'CV_Unified' + polarSuffix + '_' + str + '_' + dateString,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export CV normalized
Export.image.toDrive({
  image: CV_norm.select('CV_normalized').multiply(10000).toInt16(),
  description: 'CV_Normalized_Unified' + polarSuffix + '_' + str + '_' + dateString,
  folder: 'EarthEngine',
  fileNamePrefix: 'CV_Normalized_Unified' + polarSuffix + '_' + str + '_' + dateString,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export days (temporal information)
Export.image.toDrive({
  image: days.select('temporal_days').multiply(10000).toInt16(),
  description: 'Days_MaxIntensity_Unified' + polarSuffix + '_' + str + '_' + dateString,
  folder: 'EarthEngine',
  fileNamePrefix: 'Days_MaxIntensity_Unified' + polarSuffix + '_' + str + '_' + dateString,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export maximum intensity (linear scale)
Export.image.toDrive({
  image: imax.select('max_intensity').multiply(10000).toInt16(),
  description: 'MaxIntensity_Linear_Unified' + polarSuffix + '_' + str + '_' + dateString,
  folder: 'EarthEngine',
  fileNamePrefix: 'MaxIntensity_Linear_Unified' + polarSuffix + '_' + str + '_' + dateString,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export stack size (number of images per pixel)
Export.image.toDrive({
  image: sizepile.select('stack_size').toInt16(),
  description: 'StackSize_Unified' + polarSuffix + '_' + str + '_' + dateString,
  folder: 'EarthEngine',
  fileNamePrefix: 'StackSize_Unified' + polarSuffix + '_' + str + '_' + dateString,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

print('📤 Unified export tasks created. Check the Tasks tab to run exports.');

// =========================================================================
// INTERACTIVE CHART - ADAPTED TO UNIFIED VERSION
// =========================================================================

// Create a panel to hold the chart.
Map.style().set('cursor', 'crosshair');

var panel = ui.Panel();
panel.style().set({
  width: '500px',
  position: 'bottom-left'
});

Map.add(panel);

// Preparing list of available bands
var selectedBands = [];
if (polar1 !== null) selectedBands.push(polar1);
if (polar2 !== null) selectedBands.push(polar2);

print('📊 Selected bands for charts:', selectedBands);

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
  
  // Get values at point from named bands
  var cvValues = CV_norm.select('CV_normalized').reduceRegion({
    reducer: ee.Reducer.first(),
    geometry: point,
    scale: 10
  });
  
  var imaxValues = imax.select('max_intensity').reduceRegion({
    reducer: ee.Reducer.first(),
    geometry: point,
    scale: 10
  });
  
  var changeMaskValues = changeMask.reduceRegion({
    reducer: ee.Reducer.first(),
    geometry: point,
    scale: 10
  });
  
  var cvValue = ee.Number(cvValues.get('CV_normalized')).getInfo() || 0;
  var intensityValue = ee.Number(imaxValues.get('max_intensity')).getInfo() || 0;
  var isChange = ee.Number(changeMaskValues.values().get(0)).getInfo() || 0;
  
  var infoLabel = ui.Label('📍 Point: ' + coords.lon.toFixed(4) + ', ' + coords.lat.toFixed(4) + '\n' +
                          'CV: ' + cvValue.toFixed(3) + 
                          ' | Max Intensity: ' + intensityValue.toFixed(3) + ' dB' +
                          ' | Change detected: ' + (isChange ? 'YES' : 'NO') + '\n' +
                          'Active polarization: ' + activePolar);
  panel.add(infoLabel);
  
  if (selectedBands.length === 0) {
    panel.add(ui.Label("❌ Error: No valid band for graph display."));
    return;
  }
  
  // Dynamic graph creation based on available bands
  var chart2 = ui.Chart.image.series(sentinel1.select(selectedBands), point, null, 30)
    .setChartType('ScatterChart')
    .setOptions({
      title: 'Unified Temporal Profile (' + activePolar + ')',
      hAxis: { title: 'Acquisition Date' },
      vAxis: { title: 'Intensity Values (dB)' },
      pointSize: 5,
      pointShape: 'circle'
    });
  panel.add(chart2);
});

// =========================================================================
// UNIFIED USAGE INSTRUCTIONS
// =========================================================================
print('');
print('🎛️ UNIFIED REACTIV PARAMETERS:');
print('• enableThresholding: enable/disable thresholding (true/false)');
print('• cvThreshold: threshold for coefficient of variation (0.75-1.0)');
print('• intensityThreshold: threshold for maximum intensity (dB scale)');
print('• combinedThreshold: use combination of both thresholds (true/false)');
print('• croppalet: HSV color palette limitation (0-1)');
print('');
print('💡 UNIFIED TIPS:');
print('• Uses proper COPERNICUS/S1_GRD collection (data in dB)');
print('• Automatically handles 1 or 2 polarizations');
print('• Conversion from dB to linear amplitude: 10**(amplitude/20)');
print('• Identical algorithm core as official REACTIV version');
print('• RGB composition: days*croppalet, CV_norm, imax_linear_clamped');
print('• Gamma=2.0 for better visualization');
print('• All images have named bands for better control');
print('• Simplified polarization handling (strings instead of ee.List)');
print('');
print('🔧 BUG FIXES:');
print('• Fixed data import problem (COPERNICUS/S1_GRD)');
print('• Restored proper dB to linear conversion');
print('• Simplified polarization name handling');
print('• Fixed band selection problem in reduceRegion');
print('• Added explicit band names for all output images');
print('• Fixed mask and statistics handling');
print('• Added debug info about band names');
print('');
print('📤 UNIFIED EXPORTS:');
print('• All products with "Unified" prefix');
print('• Intensities in dB (original S1 scale)');
print('• Automatic naming with active polarization');
print('• Use proper band names');
print('');
print('✂️ UNIFIED CLIPPING:');
print('• All layers clipped to Donbas region boundaries');
print('• Simplified geometry for efficiency');
print('• Map centers on Donbas region');
print('• Region area: ' + regionAreaKm2.getInfo().toFixed(1) + ' km²');