/*
==================================================================================
📌 SEASONAL REACTIV for SAR Change Detection
==================================================================================
Description:
- This script implements the REACTIV (Radar for Evaluating Actively Cultural Areas 
  using Temporal Information and Visualization) algorithm with seasonal filtering
- Focuses analysis on specific seasonal periods (e.g., spring growing season)
  across multiple years to detect systematic changes in land use patterns
- Uses temporal coefficient of variation to identify areas of change while
  filtering out seasonal vegetation cycles and focusing on persistent changes
- Seasonal filtering helps distinguish between natural phenological cycles
  and anthropogenic changes (urban development, agricultural intensification, etc.)

Based on: Elise Colin Koeniguer et al. REACTIV methodology
References:
- Colin Koeniguer, E., et al. (2018). "Colored visualization of multitemporal SAR 
  data for change detection: issues and methods." EUSAR 2018.
- Colin Koeniguer, E., et al. (2018). "Visualisation des changements sur séries 
  temporelles radar : méthode REACTIV évaluée à l'échelle mondiale sous Google Earth Engine"
- Temporal filtering methodology for seasonal analysis in SAR time series

==================================================================================
*/

// =========================================================================
// 1. User parameters 
// =========================================================================
var str1='2022-01-01';     // Start of the Observation 'YY-MM-dd'
var str2='2025-12-31';     // End of the Observation 'YY-MM-dd'
var seasonStart = '04-01'; // Start of the season (April 1) 'MM-dd'
var seasonEnd = '06-30';   // End of the season (June 30) 'MM-dd' 

// NEW THRESHOLDING PARAMETERS
var enableThresholding = true;        // Enable/disable thresholding
var cvThreshold = 0.25;              // Threshold for coefficient of variation (0.1-0.5)
var intensityThreshold = 0.02;       // Threshold for maximum intensity
var combinedThreshold = true;        // Use combination of both thresholds (true/false)

// Check region size and complexity
var regionArea = Donbas.geometry().area();
var regionAreaKm2 = regionArea.divide(1000000); // Convert to km²
print('📐 Donbas region area (km²):', regionAreaKm2);

// If region is too large, consider simplifying or processing in tiles
var maxAreaKm2 = 50000; // Maximum recommended area
var isRegionTooLarge = regionAreaKm2.gt(maxAreaKm2);
print(⚠️ Region might be too large for processing:', isRegionTooLarge);

// Simplify geometry if needed to reduce processing complexity
var simplifiedDonbas = Donbas.map(function(feature) {
  return feature.simplify(100); // Simplify to 100m tolerance
});

// Study area: Donbas region
var geometry = simplifiedDonbas.geometry();

// Center the map on the Donbas region
Map.centerObject(simplifiedDonbas, 13);
Map.addLayer(simplifiedDonbas, {color: 'red'}, 'Donbas Region', false);

// ----------------------------------------------------------------------
// COMPUTE COLLECTION
// ----------------------------------------------------------------------

var startDate = ee.Date(str1);
var endDate = ee.Date(str2);
var ds = endDate.difference(startDate, 'day');

// Load the Sentinel-1 ImageCollection filtered by Donbas region
var sentinel1_liste = ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT')
  .filterDate(startDate, endDate)
  .filterBounds(simplifiedDonbas);

var NbOrbit = sentinel1_liste.aggregate_count_distinct('relativeOrbitNumber_start');
var ListOrbits = sentinel1_liste.aggregate_array('relativeOrbitNumber_start');
// find orbit numbers and their frequency
var freq = ee.Dictionary(ee.List(ListOrbits).reduce(ee.Reducer.frequencyHistogram()));
var array = ee.Array([freq.keys().map(ee.Number.parse), freq.values()]);
// orbit choice : first, the one with the max frequency
var frequences = array.slice(0,-1);
var arraysort = array.sort(frequences);
var index = ee.Number(NbOrbit).add(-1);
var orbite = arraysort.get([0,ee.Number(index)]);
// find images with the chosen orbit
var sentinel1 = sentinel1_liste.filterMetadata('relativeOrbitNumber_start', 'equals', orbite);

// Ascending or Descending
var modePass = sentinel1.aggregate_array('orbitProperties_pass').distinct();
// Acquisition Mode 
var acquisitionMode = sentinel1.aggregate_array('instrumentMode').distinct();

// Polarimetric channels
var allPolarizations = sentinel1.aggregate_array('transmitterReceiverPolarisation');
// Intersection of common polarimetric channels for all orbits
var commonPolarizations = allPolarizations.iterate(function(list, intersection) {
  return ee.List(intersection).removeAll(ee.List(intersection).removeAll(ee.List(list)));
}, allPolarizations.get(0));  // Initialize with the first list
var polarList = ee.List(commonPolarizations);
var NumberOfChannels = polarList.length();
var polar1 = ee.Algorithms.If(NumberOfChannels.gte(1), polarList.get(0), null);
var polar2 = ee.Algorithms.If(NumberOfChannels.eq(2), polarList.get(1), null);
var polar1 = ee.List([polar1]); 
var polar2 = ee.List([polar2]); 

// ----------------------------------------------------------------------
// COLLECTION SUMMARY
// ----------------------------------------------------------------------
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
  '📅 Full Date Range': str1 + ' to ' + str2,
  '🌸 Season': seasonStart + ' to ' + seasonEnd,
  '🎚️ CV Threshold': cvThreshold,
  '📊 Intensity Threshold': intensityThreshold,
  '🔧 Thresholding Enabled': enableThresholding,
  '👻 Transparency Mode': 'ON - non-change areas transparent'
});

// Clean display in a single print()
print('📊 Sentinel-1 Seasonal Summary:', summary);

// ----------------------------------------------------------------------
// Does the season overlap two years? 
// ----------------------------------------------------------------------

// Function to convert 'MM-dd' to day of year
function getDayOfYear(dateString, year) {
  var date = ee.Date.fromYMD(year, ee.Number.parse(dateString.split('-')[0]), ee.Number.parse(dateString.split('-')[1]));
  return date.getRelative('day', 'year');
}
// Utility function to separate MM-dd
function parseMonthDay(mmddString) {
  var parts = ee.String(mmddString).split('-');
  return {
    month: ee.Number.parse(parts.get(0)),
    day:   ee.Number.parse(parts.get(1))
  };
}
var startMD = parseMonthDay(seasonStart); // {month: x, day: y}
var endMD   = parseMonthDay(seasonEnd);   // {month: x, day: y}
// Test to see if the season "overflows" into the following year
var bridging = endMD.month.lt(startMD.month)
  .or(
    endMD.month.eq(startMD.month).and(endMD.day.lt(startMD.day))
  );

// ----------------------------------------------------------------------
// Calculation of years to cover
// ----------------------------------------------------------------------
var startYear = startDate.get('year');
var endYear   = endDate.get('year');
// List of years [startYear .. endYear]
var years = ee.List.sequence(startYear, endYear);
// Creation of the season interval for a year y 
function makeSeasonInterval(y) {
  y = ee.Number(y);
  // If bridging is true, the end of season is year y+1, otherwise it's the same year y  
  var sYear = y;
  var eYear = ee.Algorithms.If(bridging.eq(1), y.add(1), y);
  // Build the startSeason and endSeason for this year
  var startSeason = ee.Date.fromYMD(sYear, startMD.month, startMD.day);
  var endSeason   = ee.Date.fromYMD(eYear, endMD.month, endMD.day);
  // Calculate overlap with [startDate, endDate] 
  // => take max for start, min for end, in terms of timestamps
  var overlapStart = ee.Date(
    ee.Number(startSeason.millis()).max(startDate.millis())
  );
  var overlapEnd = ee.Date(
    ee.Number(endSeason.millis()).min(endDate.millis())
  );
  // Returns a dictionary describing the interval
  return ee.Dictionary({
    'year':      y,
    'start': overlapStart.format("YYYY-MM-dd'T'HH:mm:ss'Z'"),
    'end':   overlapEnd.format("YYYY-MM-dd'T'HH:mm:ss'Z'"),
  });
}

// Generate list of season intervals for each year
var rawIntervals = years.map(function(y) {
  return makeSeasonInterval(y);
});

var rawClient = rawIntervals.getInfo(); // => Array of dictionary objects
var filteredClient = rawClient.filter(function(d) {
  // d.startDate and d.endDate are strings "YYYY-MM-dd..."
  var startMs = new Date(d.start).getTime();
  var endMs   = new Date(d.end).getTime();
  return startMs < endMs;  // local boolean
});
var validIntervals = ee.List(filteredClient.map(function(obj) {
  return ee.Dictionary(obj);
}));

/****  Add fractions [0..1] to each interval ****/
var finalIntervals = validIntervals.map(function(d) {
  d = ee.Dictionary(d);
  var st = ee.Date(d.get('start'));
  var en = ee.Date(d.get('end'));
  var startFrac = st.difference(startDate, 'day').divide(ds); // startFrac = (st - startDate) / ds
  var endFrac   = en.difference(startDate, 'day').divide(ds); // endFrac   = (en - startDate) / ds
  // Bound within [0,1], just in case
  var sf = startFrac.max(0).min(1);
  var ef = endFrac.max(0).min(1);
  return d.combine({
    'startFrac': sf,
    'endFrac':   ef
  });
});

// ----------------------------------------------------------------------
// SEASONAL FILTERING
// ----------------------------------------------------------------------

function filterByFinalIntervals(imageCollection, finalIntervals) {
    return ee.ImageCollection(
        finalIntervals.iterate(function(interval, prev) {
            interval = ee.Dictionary(interval);
            var start = ee.Date(interval.get('start'));
            var end = ee.Date(interval.get('end'));
            var filtered = imageCollection.filterDate(start, end);
            return ee.ImageCollection(prev).merge(filtered);
        }, ee.ImageCollection([]))
    );
}

// Apply the filter function
var filteredSentinel1 = filterByFinalIntervals(sentinel1, finalIntervals);
// Result display
print('📂 Sentinel-1 filtered collection:', filteredSentinel1);
sentinel1=filteredSentinel1;

// ----------------------------------------------------------------------
// DATES 
// ----------------------------------------------------------------------
var range = sentinel1.select(polar1).reduceColumns(ee.Reducer.minMax(), ["system:time_start"]);
var minDate = ee.Date(range.get('min')).format('YYYY-MM-dd');
var maxDate = ee.Date(range.get('max')).format('YYYY-MM-dd');
// Format dates to make them readable
var okMap2 = sentinel1.select(polar1).map(function(image) {
  return image.set('date', image.date().format('YYYY-MM-dd'));
});
// Obtain a single list of sorted dates
var datesList = okMap2.aggregate_array('date').distinct().sort();
// Count the number of unique dates
var datesCount = datesList.size();
// Formatted summary of dates
var dateSummary = ee.Dictionary({
  '📅 Date Range': minDate.cat(' → ').cat(maxDate),
  '📊 Number of Unique Dates': datesCount,
  '🗓️ Number of Seasons': validIntervals.size()
});
// Compact display
print('📆 Sentinel-1 Date Summary:', dateSummary);

// ----------------------------------------------------------------------
// REACTIV RGB COMPUTATION
// ----------------------------------------------------------------------

// This function applies to each image the linear scale
var amplitude = function(image) {
  var imlin = image.expression(
    'sqrt(intensity)', {
      'intensity': image
  });
  return imlin; // conversion in linear, then compute mean: classical mean
};

// **2 polarizations**
var result = ee.Algorithms.If(
  NumberOfChannels.eq(2),
  
  // If NumberOfChannels == 2
  (function() {
    var stdLinear1 = sentinel1.select(polar1).map(amplitude).reduce(ee.Reducer.stdDev());
    var meanLinear1 = sentinel1.select(polar1).map(amplitude).reduce(ee.Reducer.mean());
    var magic1 = stdLinear1.divide(meanLinear1);

    var stdLinear2 = sentinel1.select(polar2).map(amplitude).reduce(ee.Reducer.stdDev());
    var meanLinear2 = sentinel1.select(polar2).map(amplitude).reduce(ee.Reducer.mean());
    var magic2 = stdLinear2.divide(meanLinear2);
    
    var imax1 = sentinel1.select(polar1).max();
    var imax2 = sentinel1.select(polar2).max();
    var magic = magic1.max(magic2);
    var imax = imax1.max(imax2);

    // Time calculation function for first polarization
    var time1 = function(image) {
      var days = image.date().difference(startDate, 'day').divide(ds);
      return image.where(image.lt(imax1), 0).where(image.gte(imax1), days);
    };
    var days1 = sentinel1.select(polar1).map(time1).sum();

    // Time calculation function for second polarization
    var time2 = function(image) {
      var days = image.date().difference(startDate, 'day').divide(ds);
      return image.where(image.lt(imax2), 0).where(image.gte(imax2), days);
    };
    var days2 = sentinel1.select(polar2).map(time2).sum();

    var days = days2.where(magic2.lte(magic1), days1);
    
    return {
      days: days,
      magic: magic,
      imax: imax,
      magic1: magic1,
      magic2: magic2
    };
  })(),
  
  // If NumberOfChannels == 1
  (function() {
    var stdLinear1 = sentinel1.select(polar1).map(amplitude).reduce(ee.Reducer.stdDev());
    var meanLinear1 = sentinel1.select(polar1).map(amplitude).reduce(ee.Reducer.mean());
    var magic1 = stdLinear1.divide(meanLinear1);
    
    var imax1 = sentinel1.select(polar1).max();
    var magic = magic1;
    var imax = imax1;

    // Time calculation function for first and single polarization
    var time1 = function(image) {
      var days = image.date().difference(startDate, 'day').divide(ds);
      return image.where(image.lt(imax1), 0).where(image.gte(imax1), days);
    };
    var days1 = sentinel1.select(polar1).map(time1).sum();

    return {
      days: days1,
      magic: magic,
      imax: imax,
      magic1: magic1,
      magic2: null
    };
  })()
);

// Extraction of final Results
var days = ee.Image(ee.Dictionary(result).get('days'));
var magic = ee.Image(ee.Dictionary(result).get('magic'));
var imax = ee.Image(ee.Dictionary(result).get('imax'));
var magic1 = ee.Image(ee.Dictionary(result).get('magic1'));
var magic2 = ee.Dictionary(result).get('magic2');

// Clip all results to Donbas region
days = days.clip(simplifiedDonbas);
magic = magic.clip(simplifiedDonbas);
imax = imax.clip(simplifiedDonbas);
magic1 = magic1.clip(simplifiedDonbas);
if (magic2 !== null) {
  magic2 = ee.Image(magic2).clip(simplifiedDonbas);
}

// Images of Number of images: StackSize
var unit = function(image) {
  var imunit = image.multiply(0).add(1);
  return imunit; // conversion in linear, then compute mean: classical mean
};
var sizepile=sentinel1.select(polar1).map(unit).sum().clip(simplifiedDonbas); 

// Parameter for dynamics
var mu=0.2286; // Theoretical mean for Rayleigh Nakagami L=4.9
var stdmu=ee.Image(0.1616);
var stdmu=stdmu.divide(sizepile.sqrt()); // Theoretical std for Rayleigh Nakagami L=4.9
var magicnorm=magic.subtract(mu).divide(stdmu.multiply(10)).clamp(0,1);

var rgb=ee.Image.cat(days,magicnorm,(imax).multiply(1).clamp(0,1)).hsvToRgb().clip(simplifiedDonbas);

// =========================================================================
// SEASONAL THRESHOLDING - NEW FUNCTIONALITY
// =========================================================================

// Creating thresholding masks
var cvMask = magicnorm.gte(cvThreshold);
var intensityMask = imax.gte(intensityThreshold);

// Combining masks depending on settings
var changeMask;
if (combinedThreshold) {
  changeMask = cvMask.and(intensityMask);
  print('🎭 Using combined CV and intensity thresholding for seasonal analysis');
} else {
  changeMask = cvMask;
  print('🎭 Using CV thresholding only for seasonal analysis');
}

// Creating visualization with thresholding
var rgbThresholded;
if (enableThresholding) {
  // Apply transparency mask - only changes will be visible
  rgbThresholded = rgb.updateMask(changeMask);
  
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
  
  print('📈 Seasonal Thresholding Statistics:');
  print('Total pixels:', totalPixels);
  print('Change pixels:', changePixels);
  print('🔍 Areas without changes will be transparent');
  
} else {
  rgbThresholded = rgb;
  print('🚫 Thresholding disabled - showing all seasonal data');
}

// =========================================================================
// VISUALIZATION
// =========================================================================

Map.setOptions('satellite');
var visparams = {min: [0, 0, 0],
                 max: [1, 1, 1],
                 gamma: 2};

// Main layers
Map.addLayer(rgb, visparams,'REACTIV Seasonal Original', false);
Map.addLayer(rgbThresholded, visparams,'REACTIV Seasonal Transparent');

// Add additional layers
Map.addLayer(magic1.clamp(0,0.5).divide(0.5).pow(2), {min: 0, max: 1, gamma: 1}, 'CV, Polar 1', false);
if (magic2 !== null) {
  Map.addLayer(ee.Image(magic2).clamp(0,0.5).divide(0.5).pow(2), {min: 0, max: 1, gamma: 1}, 'CV, Polar 2', false);
}

// Add masks as additional layers for control
Map.addLayer(changeMask.selfMask(), {palette: ['white']}, 'Seasonal Change Mask', false);
Map.addLayer(cvMask.selfMask(), {palette: ['yellow']}, 'CV Mask', false);
Map.addLayer(intensityMask.selfMask(), {palette: ['cyan']}, 'Intensity Mask', false);

// =========================================================================
// Export Options
// =========================================================================

// Get the date string for file naming
var dateString = startDate.format('YYYY-MM-dd').getInfo() + '_to_' + endDate.format('YYYY-MM-dd').getInfo();
var seasonString = seasonStart.replace('-', '') + '_to_' + seasonEnd.replace('-', '');
var thresholdSuffix = enableThresholding ? '_thresh_' + cvThreshold.toString().replace('.', '') : '';

// Export REACTIV RGB visualization (original)
Export.image.toDrive({
  image: rgb.multiply(255).toByte(),
  description: 'REACTIV_Seasonal_RGB_' + seasonString + '_' + dateString,
  folder: 'GEE_exports',
  fileNamePrefix: 'REACTIV_Seasonal_RGB_' + seasonString + '_' + dateString,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export REACTIV RGB visualization (with thresholding and transparency)
Export.image.toDrive({
  image: rgbThresholded.multiply(255).toByte(),
  description: 'REACTIV_Seasonal_RGB_Transparent_' + seasonString + '_' + dateString + thresholdSuffix,
  folder: 'GEE_exports',
  fileNamePrefix: 'REACTIV_Seasonal_RGB_Transparent_' + seasonString + '_' + dateString + thresholdSuffix,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export change mask
Export.image.toDrive({
  image: changeMask.toByte(),
  description: 'Seasonal_Change_Mask_' + seasonString + '_' + dateString + thresholdSuffix,
  folder: 'GEE_exports',
  fileNamePrefix: 'Seasonal_Change_Mask_' + seasonString + '_' + dateString + thresholdSuffix,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export CV Polar 1
Export.image.toDrive({
  image: magic1.multiply(10000).toInt16(),
  description: 'CV_Polar1_Seasonal_' + seasonString + '_' + dateString,
  folder: 'GEE_exports',
  fileNamePrefix: 'CV_Polar1_Seasonal_' + seasonString + '_' + dateString,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export CV Polar 2 (if exists)
if (magic2 !== null) {
  Export.image.toDrive({
    image: ee.Image(magic2).multiply(10000).toInt16(),
    description: 'CV_Polar2_Seasonal_' + seasonString + '_' + dateString,
    folder: 'GEE_exports',
    fileNamePrefix: 'CV_Polar2_Seasonal_' + seasonString + '_' + dateString,
    region: simplifiedDonbas.geometry(),
    scale: 10,
    crs: 'EPSG:4326',
    maxPixels: 1e13
  });
}

// Export days (temporal information)
Export.image.toDrive({
  image: days.multiply(10000).toInt16(),
  description: 'Days_MaxIntensity_Seasonal_' + seasonString + '_' + dateString,
  folder: 'GEE_exports',
  fileNamePrefix: 'Days_MaxIntensity_Seasonal_' + seasonString + '_' + dateString,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export maximum intensity
Export.image.toDrive({
  image: imax.multiply(10000).toInt16(),
  description: 'MaxIntensity_Seasonal_' + seasonString + '_' + dateString,
  folder: 'GEE_exports',
  fileNamePrefix: 'MaxIntensity_Seasonal_' + seasonString + '_' + dateString,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export normalized CV (magicnorm)
Export.image.toDrive({
  image: magicnorm.multiply(10000).toInt16(),
  description: 'CV_Normalized_Seasonal_' + seasonString + '_' + dateString,
  folder: 'GEE_exports',
  fileNamePrefix: 'CV_Normalized_Seasonal_' + seasonString + '_' + dateString,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export stack size (number of images per pixel)
Export.image.toDrive({
  image: sizepile.toInt16(),
  description: 'StackSize_Seasonal_' + seasonString + '_' + dateString,
  folder: 'GEE_exports',
  fileNamePrefix: 'StackSize_Seasonal_' + seasonString + '_' + dateString,
  region: simplifiedDonbas.geometry(),
  scale: 10,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

print('📤 Export tasks created. Check the Tasks tab to run exports.');

// ----------------------------------------------------------------------
// COLORED TIMELINE LEGEND
// ----------------------------------------------------------------------

var clientList = finalIntervals.getInfo();

function fractionToHex(hueFrac) {
  hueFrac = ee.Number(hueFrac).clamp(0, 1);
  var hsvImg = ee.Image.cat([hueFrac, 1, 1]).rename(['h', 's', 'v']);
  var rgbImg = hsvImg.hsvToRgb().multiply(255).round().toInt();
  var stats = rgbImg.reduceRegion({
    reducer: ee.Reducer.first(),
    geometry: ee.Geometry.Point([0, 0]), // Arbitrary point
    scale: 10,
    bestEffort: true,
    maxPixels: 1e6
  });
  var r = ee.Number(stats.get('red')).format('%02X');
  var g = ee.Number(stats.get('green')).format('%02X');
  var b = ee.Number(stats.get('blue')).format('%02X');
  return ee.String('#').cat(r).cat(g).cat(b);
}

var withHex = finalIntervals.map(function(elem) {
  elem = ee.Dictionary(elem);
  var sf = ee.Number(elem.get('startFrac'));
  var ef = ee.Number(elem.get('endFrac'));
  var startHex = fractionToHex(sf);
  var endHex   = fractionToHex(ef);
  return elem.combine({
    'startHex': startHex,
    'endHex':   endHex
  });
});
var clientList = withHex.getInfo();

function makeLegend(vis,d1,d2) {
  var lon = ee.Image.pixelLonLat().select('longitude');
  var gradient = lon.multiply(1/360.0);
  var legendImage = gradient.visualize(vis);
  var thumb = ui.Thumbnail({
    image: legendImage,
    params: {bbox:'0,0,360,4', dimensions:'70x10'},
    style: {padding: '1px', position: 'bottom-center',backgroundColor:'white'}
    });
  var panel = ui.Panel({
  widgets: [
  ui.Label(d1),
  ui.Label(d2),
  ],
  layout: ui.Panel.Layout.flow('vertical'),
  style: {stretch: 'horizontal',backgroundColor:'white',color:'blue'}
  });
  return ui.Panel({style: {backgroundColor: 'white'}}).add(panel).add(thumb);
}

function formatDate(dateString) {
  return dateString.slice(2, 10); // Keep only "YY-MM-dd".
}

for (var i=0; i<clientList.length; i++) {
  var startHex = clientList[i].startHex;
  var endHex   = clientList[i].endHex;
  var numElements = 10; // Set the number of elements in the palette - Adjust the number to suit your needs
  
  var startFrac = clientList[i].startFrac;
  var endFrac = clientList[i].endFrac;

  // Automatic step calculation based on number of elements
  var step = (endFrac - startFrac) / (numElements - 1);
  
  // Fraction list creation
  var fracList = ee.List.sequence(startFrac, endFrac, step);

  // Generate palette by applying fractionToHex to each fracList element
  var palettehsv = fracList.map(function(hueFrac) {
    return ee.String(fractionToHex(hueFrac));
  });

  var pal = {
    min: 0,
    max: 1,
    palette: palettehsv
  };
    // Format dates before passing them to makeLegend
  var formattedStartDate = formatDate(clientList[i].start);
  var formattedEndDate = formatDate(clientList[i].end);
  Map.add(makeLegend(pal, formattedStartDate, formattedEndDate));
}

// ----------------------------------------------------------------------
// TIME PROFILE PLOT
// ----------------------------------------------------------------------
// Create the title label.
var title = ui.Label('Click to inspect (Seasonal Thresholding: ' + (enableThresholding ? 'ON' : 'OFF') + ')');
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
  var pointValues = ee.Image.cat(magicnorm, imax, changeMask).reduceRegion({
    reducer: ee.Reducer.first(),
    geometry: point,
    scale: 10
  }).getInfo();
  
  var cvValue = pointValues.cv || 0;
  var intensityValue = pointValues.intensity || 0;
  var isChange = pointValues.cv_1 || 0;
  
  var infoLabel = ui.Label('CV: ' + cvValue.toFixed(3) + 
                          ' | Intensity: ' + intensityValue.toFixed(3) + 
                          ' | Seasonal change detected: ' + (isChange ? 'YES' : 'NO'));
  panel.add(infoLabel);
  
  print('Coordinates for the plot: ', point);
  if (selectedBandsJS.length === 0) {
    if (panel.widgets().length() === 0) { // Prevent multiple display of message
      panel.add(ui.Label("❌ Error: No valid band for graph display."));
    }
    return; // Stop execution here
  }
  // Dynamic graph creation based on available bands
  var chart2 = ui.Chart.image.series(sentinel1.select(selectedBandsJS), point, null, 30)
    .setChartType('ScatterChart')
    .setOptions({
      title: 'Seasonal Temporal Profile at: ' + coords.lon.toFixed(4) + ', ' + coords.lat.toFixed(4),
      hAxis: { title: 'Acquisition Date' },
      vAxis: { title: 'Intensity Values (linear)' },
    });

  panel.add(chart2);
});

// =========================================================================
// SEASONAL USAGE INSTRUCTIONS
// =========================================================================
print('');
print('🌸 SEASONAL THRESHOLDING PARAMETERS:');
print('• enableThresholding: enable/disable thresholding (true/false)');
print('• cvThreshold: threshold for coefficient of variation (0.1-0.5)');
print('• intensityThreshold: threshold for maximum intensity');
print('• combinedThreshold: use combination of both thresholds (true/false)');
print('• seasonStart/seasonEnd: seasonal analysis period');
print('');
print('💡 SEASONAL TIPS:');
print('• Thresholding on seasonal data shows changes in selected period');
print('• Increase cvThreshold to see only stronger seasonal changes');
print('• Decrease cvThreshold to see more subtle changes');
print('• Experiment with combinedThreshold for different strategies');
print('• Check auxiliary masks in layer panel for debugging');
print('• Areas without changes will be transparent - underlying layers will be visible');
print('• Color legends show temporal distribution within selected season');
print('');
print('🔬 SEASONAL METHODOLOGY:');
print('• Focuses analysis on specific seasonal windows across multiple years');
print('• Helps distinguish natural phenological cycles from anthropogenic changes');
print('• Useful for agricultural monitoring, construction seasons, conflict periods');
print('• Cross-year seasonal comparison reveals systematic trends');
print('• HSV visualization: Hue=timing, Saturation=variability, Value=intensity');
print('');
print('📊 INTERPRETATION:');
print('• Consistent seasonal changes indicate systematic land use shifts');
print('• High CV in seasonal windows suggests disturbance or development');
print('• Timing information shows when within season changes occurred');
print('• Multi-year patterns reveal long-term trends vs. episodic events');
print('');
print('✅ Seasonal REACTIV processing complete!');