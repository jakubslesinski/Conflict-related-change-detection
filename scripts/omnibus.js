/*
==================================================================================
📌 SEQUENTIAL OMNIBUS ALGORITHM for SAR Change Detection
==================================================================================
Description:
- This script implements the Omnibus test for change detection in a time sequence 
  of polarimetric SAR data
- The omnibus test is a multivariate statistical test that can detect changes 
  in the mean and covariance structure of polarimetric SAR data
- Sequential processing allows identification of multiple change points in time series
- The algorithm outputs: total number of changes (cmap), first change (smap), 
  last change (fmap), and binary change maps for each date (bmap)

Author: Based on Morton Canty's GEE implementation
Reference:
A. A. Nielsen, K. Conradsen and H. Skriver, "Omnibus test for change detection in a time sequence of polarimetric SAR data," 2016 IEEE International Geoscience and Remote Sensing Symposium (IGARSS), Beijing, China, 2016, pp. 3398-3401, doi: 10.1109/IGARSS.2016.7729878.

Additional References:
- Conradsen, K., Nielsen, A. A., Schou, J., & Skriver, H. (2003). A test statistic in the complex Wishart distribution and its application to change detection in polarimetric SAR data. IEEE Transactions on Geoscience and Remote Sensing, 41(1), 4-19.
- Morton J. Canty, "Image Analysis, Classification and Change Detection in Remote Sensing", 4th Edition, CRC Press 2019

==================================================================================
*/

// =========================================================================
// Sequential Omnibus Algorithm - SIMPLIFIED VERSION
// Modified for Donbas Region with Export Options
// 
// SIMPLIFIED VERSION - No seasonal filtering, just standard time series analysis
// =========================================================================

// -------------------------------------------------------------------------
// LOAD DONBAS REGION
// -------------------------------------------------------------------------
// IMPORTANT: You need to either:
// 1. Import your Donbas shapefile as an asset and define it here:
//    var Donbas = ee.FeatureCollection('users/yourUsername/Donbas');
// 2. OR have it already loaded in your workspace before running this script
// -------------------------------------------------------------------------

// Check if Donbas is defined, if not provide instructions
if (typeof Donbas === 'undefined') {
  print('⚠️ ERROR: Donbas region not found!');
  print('Please either:');
  print('1. Import your Donbas shapefile asset at the top of this script:');
  print('   var Donbas = ee.FeatureCollection("users/yourUsername/Donbas");');
  print('2. OR load it in your workspace before running this script');
  throw new Error('Donbas region not defined');
}

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

// Define styling for the shapefile  
var shp_styl = {
  color: 'yellow',
  fillColor: '00000000',
};

Map.addLayer(simplifiedDonbas.style(shp_styl), {}, 'Donbas Region');

// Use simplified Donbas geometry
var geometry = simplifiedDonbas.geometry();

Map.setOptions('satellite');
Map.style().set('cursor', 'crosshair');
Map.centerObject(simplifiedDonbas, 8);

// ******************************************
// Front End for Sequential Omnibus Algorithm
// ******************************************
var omb = require('users/mortcanty/changedetection:omnibus');
var util = require('users/mortcanty/changedetection:utilities');
var jet = ['black','blue','cyan', 'yellow','red'];

// =========================================================================
// 1. User Parameters - SIMPLIFIED
// =========================================================================
var significance = 0.0001;   
var relorbitnumber = 'any'; 
var orbitpass = 'ASCENDING';
var median = true;

// Date range
var startDate = '2022-01-01';  
var endDate = '2022-12-31';    

// Asset export path
var shortDateRange = ee.Date(startDate).format('yyMMdd').getInfo() + '_' + ee.Date(endDate).format('yyMMdd').getInfo();

// =========================================================================
// 2. Collection Filtering
// =========================================================================
var bnds = ee.List(geometry.bounds().coordinates().get(0));
var p0 = ee.Geometry.Point(bnds.get(0));
var p1 = ee.Geometry.Point(bnds.get(1));
var p2 = ee.Geometry.Point(bnds.get(2));
var p3 = ee.Geometry.Point(bnds.get(3));

var collection = ee.ImageCollection('COPERNICUS/S1_GRD') 
                   .filterBounds(p0) 
                   .filterBounds(p1)
                   .filterBounds(p2) 
                   .filterBounds(p3) 
                   .filterDate(ee.Date(startDate), ee.Date(endDate)) 
                   .filter(ee.Filter.eq('transmitterReceiverPolarisation', ['VV','VH'])) 
                   .filter(ee.Filter.eq('resolution_meters', 10)) 
                   .filter(ee.Filter.eq('instrumentMode', 'IW'))
                   .filter(ee.Filter.eq('orbitProperties_pass', orbitpass)); 

if (relorbitnumber != 'any'){
  collection = collection.filter(ee.Filter.eq('relativeOrbitNumber_start', relorbitnumber));                        
} 

collection = collection.sort('system:time_start');  

// =========================================================================
// 3. Process Collection
// =========================================================================
var acquisition_times = ee.List(collection.aggregate_array('system:time_start'));
var count = acquisition_times.length().getInfo();

if (count === 0){ 
  print('❌ No images found');
  print('   Check your date range and orbit settings');
} else if (count > 100) {
  print('❌ ERROR: Time series too long (' + count + ' images)');
  print('   Maximum recommended: 100 images');
  print('   Current settings:');
  print('   - Date range:', startDate, 'to', endDate);
  print('   - Orbit:', relorbitnumber);
  print('   - Pass:', orbitpass);
  print('   Suggestions:');
  print('   - Use shorter date range');
  print('   - Set specific orbit number instead of "any"');
  
  // Show orbit frequency to help user choose
  var allOrbits = collection.aggregate_array('relativeOrbitNumber_start');
  var orbitFreq = ee.Dictionary(allOrbits.reduce(ee.Reducer.frequencyHistogram()));
  print('\n📊 Orbit frequency in your data:');
  print(orbitFreq);
  print('💡 TIP: Use the most frequent orbit number to reduce images.');
} else {
  print('📊 Omnibus Analysis Summary:');
  print('📅 Date Range:', startDate, 'to', endDate);
  print('🛰️ Orbit Pass:', orbitpass);
  print('🔢 Orbit Number:', relorbitnumber);
  print('📸 Number of Images:', count);
  print('📈 Significance Level:', significance);
  print('🏭 Median Filter:', median);
  
  // Warning for long time series
  if (count > 50) {
    print('⚠️ WARNING: Long time series detected (' + count + ' images)');
    print('   Processing may take several minutes.');
  }
  
  print('📅 Timestamps:');
  var timestamps = acquisition_times.getInfo().map(function(d){ return new Date(d); });
  print(timestamps);
  print('🛰️ Relative Orbit Numbers:');
  print(ee.List(collection.aggregate_array('relativeOrbitNumber_start'))); 
  
  // =========================================================================
  // 4. Create Legend
  // =========================================================================
  var vis = {min: 0, max: count, palette: jet};
  Map.add(util.makeLegend(vis));
  
  // Create a list of clipped images
  var pList = collection.map(omb.get_vvvh).toList(count);
  var first = ee.Dictionary({imlist:ee.List([]),geom:geometry});
  var imList = ee.List(ee.Dictionary(pList.iterate(omb.clipList,first)).get('imlist')); 
  
  // =========================================================================
  // 5. Run the Omnibus Algorithm
  // =========================================================================
  print('🔄 Running Omnibus algorithm...');
  var result = ee.Dictionary(omb.omnibus(imList,significance,median));
  
  // Get change maps 
  var cmap = ee.Image(result.get('cmap')).byte();
  var smap = ee.Image(result.get('smap')).byte();
  var fmap = ee.Image(result.get('fmap')).byte();
  var bmap = ee.Image(result.get('bmap')).byte();
  
  var cnames = ['cmap','smap','fmap'];
  for (var i = 1; i < count; i++){
    if (i < 10) {var label = 'bmap0';} else {var label= 'bmap';}
    cnames = cnames.concat([label.concat(i.toString())]);
  }
  cnames = cnames.concat(['background']);
  
  // Background image for video export
  var background = collection.mean()
                            .select(0)
                            .multiply(ee.Image.constant(Math.log(10.0)/10.0)).exp();
  background = background.where(background.gte(1),1).clip(geometry);
  
  // Concatenate change maps
  var cmaps = ee.Image.cat(cmap,smap,fmap,bmap,background).rename(cnames);
  
  print('✅ Omnibus algorithm completed successfully');
  
  // =========================================================================
  // 6. Export Options
  // =========================================================================
  
  print('📤 Creating export tasks...');
  
  // Get the date string for file naming - shortened format
  var dateString = ee.Date(startDate).format('yyMMdd').getInfo() + '_' + ee.Date(endDate).format('yyMMdd').getInfo();
  
  // Export to Asset (original functionality)
  //var exportTask = Export.image.toAsset({
  //  image: cmaps,
  //  assetId: assetExportId,
  //  scale: 10,
  //  maxPixels: 1e9
  //});
  
  // Export cmap (number of changes)
  Export.image.toDrive({
    image: cmap,
    description: 'Omnibus_cmap_' + orbitpass + '_' + dateString,
    folder: 'EarthEngine_Omnibus',
    fileNamePrefix: 'Omnibus_cmap_' + orbitpass + '_' + dateString,
    region: simplifiedDonbas.geometry(),
    scale: 10,
    crs: 'EPSG:4326',
    maxPixels: 1e13
  });
  
  // Export smap (first change)
  Export.image.toDrive({
    image: smap,
    description: 'Omnibus_smap_' + orbitpass + '_' + dateString,
    folder: 'EarthEngine_Omnibus',
    fileNamePrefix: 'Omnibus_smap_' + orbitpass + '_' + dateString,
    region: simplifiedDonbas.geometry(),
    scale: 10,
    crs: 'EPSG:4326',
    maxPixels: 1e13
  });
  
  // Export fmap (last change)
  Export.image.toDrive({
    image: fmap,
    description: 'Omnibus_fmap_' + orbitpass + '_' + dateString,
    folder: 'EarthEngine_Omnibus',
    fileNamePrefix: 'Omnibus_fmap_' + orbitpass + '_' + dateString,
    region: simplifiedDonbas.geometry(),
    scale: 10,
    crs: 'EPSG:4326',
    maxPixels: 1e13
  });
  
  // Export all binary change maps (individual dates)
  Export.image.toDrive({
    image: bmap,
    description: 'Omnibus_bmap_' + orbitpass + '_' + dateString,
    folder: 'EarthEngine_Omnibus',
    fileNamePrefix: 'Omnibus_bmap_' + orbitpass + '_' + dateString,
    region: simplifiedDonbas.geometry(),
    scale: 10,
    crs: 'EPSG:4326',
    maxPixels: 1e13
  });
  
  // Export all change maps combined
  Export.image.toDrive({
    image: cmaps.select(['cmap','smap','fmap']),
    description: 'Omnibus_all_changes_' + orbitpass + '_' + dateString,
    folder: 'EarthEngine_Omnibus',
    fileNamePrefix: 'Omnibus_all_changes_' + orbitpass + '_' + dateString,
    region: simplifiedDonbas.geometry(),
    scale: 10,
    crs: 'EPSG:4326',
    maxPixels: 1e13
  });
  
  // Export background (mean amplitude)
  Export.image.toDrive({
    image: background.multiply(10000).toInt16(),
    description: 'Omnibus_background_' + orbitpass + '_' + dateString,
    folder: 'EarthEngine_Omnibus',
    fileNamePrefix: 'Omnibus_background_' + orbitpass + '_' + dateString,
    region: simplifiedDonbas.geometry(),
    scale: 10,
    crs: 'EPSG:4326',
    maxPixels: 1e13
  });
  
  // Export complete analysis package
  Export.image.toDrive({
    image: cmaps,
    description: 'Omnibus_complete_' + orbitpass + '_' + dateString,
    folder: 'EarthEngine_Omnibus',
    fileNamePrefix: 'Omnibus_complete_' + orbitpass + '_' + dateString,
    region: simplifiedDonbas.geometry(),
    scale: 10,
    crs: 'EPSG:4326',
    maxPixels: 1e13
  });
  
  print('📤 Export tasks created. Check the Tasks tab to run exports.');
  
  // =========================================================================
  // 7. Display Results
  // =========================================================================
  Map.centerObject(simplifiedDonbas, 10);
  
  // Create visualizations
  var cmapVis = {min: 0, max: count, palette: jet};
  var smapVis = {min: 0, max: count, palette: jet};
  var fmapVis = {min: 0, max: count, palette: jet};
  
  Map.addLayer(cmaps.select('cmap'), cmapVis, 'cmap - Total Changes', false, 0.6);
  Map.addLayer(cmaps.select('smap'), smapVis, 'smap - First Change', true, 0.6);
  Map.addLayer(cmaps.select('fmap'), fmapVis, 'fmap - Last Change', false, 0.6);
  Map.addLayer(background, {min:0, max:1}, 'Background', false);
  
  // =========================================================================
  // 8. Interactive Information Panel
  // =========================================================================
  var panel = ui.Panel({
    style: {
      width: '350px',
      position: 'bottom-right'
    }
  });
  
  var intro = ui.Label({
    value: '📊 Omnibus Results',
    style: {fontWeight: 'bold', fontSize: '16px'}
  });
  panel.add(intro);
  
  var summaryInfo = ui.Label({
    value: '📅 ' + startDate + ' to ' + endDate + '\n' +
           '📸 ' + count + ' images analyzed\n' +
           '🛰️ ' + orbitpass + ' orbit',
    style: {fontSize: '11px', color: 'gray'}
  });
  panel.add(summaryInfo);
  
  var info = ui.Label({
    value: 'Click on the map to see change information',
    style: {fontSize: '12px'}
  });
  panel.add(info);
  
  Map.add(panel);
  
  // Add click handler
  Map.onClick(function(coords) {
    var point = ee.Geometry.Point(coords.lon, coords.lat);
    
    // Check if the clicked point is within Donbas region
    var withinDonbas = simplifiedDonbas.geometry().contains(point);
    var isInside = withinDonbas.getInfo();
    
    if (!isInside) {
      info.setValue('⚠️ Click within the Donbas region');
    } else {
      // Sample the change maps at the clicked point
      var sample = cmaps.sample({
        region: point,
        scale: 10,
        geometries: true
      });
      
      var values = sample.first();
      values.evaluate(function(result) {
        if (result) {
          var props = result.properties;
          var changeInfo = '📍 Location: ' + coords.lon.toFixed(4) + ', ' + coords.lat.toFixed(4) + '\n' +
                          '🔄 Total changes: ' + props.cmap + '\n' +
                          '🟢 First change at: ' + (props.smap > 0 ? 'Image ' + props.smap : 'No change') + '\n' +
                          '🔴 Last change at: ' + (props.fmap > 0 ? 'Image ' + props.fmap : 'No change');
          info.setValue(changeInfo);
          
          // Get date of first change if exists
          if (props.smap > 0 && props.smap <= timestamps.length) {
            var firstChangeDate = timestamps[props.smap - 1];
            info.setValue(info.getValue() + '\n📅 First: ' + firstChangeDate.toLocaleDateString());
          }
          if (props.fmap > 0 && props.fmap <= timestamps.length) {
            var lastChangeDate = timestamps[props.fmap - 1];
            info.setValue(info.getValue() + '\n📅 Last: ' + lastChangeDate.toLocaleDateString());
          }
          
          // Background value
          if (props.background) {
            info.setValue(info.getValue() + '\n🏔️ Background: ' + props.background.toFixed(3));
          }
        }
      });
    }
  });
  
  // =========================================================================
  // 9. Summary Information
  // =========================================================================
  
  print('');
  print('📊 ANALYSIS RESULTS SUMMARY:');
  print('✅ Omnibus algorithm completed successfully');
  print('📅 Period analyzed:', startDate, 'to', endDate);
  print('📸 Total images processed:', count);
  print('🎯 Significance level:', significance);
  print('📈 Expected false positive rate:', (significance * 100).toFixed(4) + '%');
  print('');
  print('🎨 VISUALIZATIONS:');
  print('• cmap (Total Changes): Shows number of changes per pixel');
  print('• smap (First Change): Shows when first change occurred');  
  print('• fmap (Last Change): Shows when last change occurred');
  print('• Background: Mean amplitude of time series');
  print('');
  print('📤 EXPORTS AVAILABLE:');
  print('• Omnibus_cmap: Number of changes per pixel');
  print('• Omnibus_smap: First change detection');
  print('• Omnibus_fmap: Last change detection');
  print('• Omnibus_bmap: Individual change dates');
  print('• Omnibus_background: Reference background');
  print('• Omnibus_all_changes: Combined change maps');
  print('• Omnibus_complete: Full analysis package');
  print('');
  print('📁 Export folder: EarthEngine_Omnibus');
  print('📏 Region: Donbas (' + regionAreaKm2.getInfo().toFixed(1) + ' km²)');
  print('');
  print('💡 INTERPRETATION TIPS:');
  print('• Higher cmap values = more dynamic areas');
  print('• smap timing = early changes (conflicts, disasters)');
  print('• fmap timing = recent changes (reconstruction, new activities)');
  print('• Compare with background to understand change magnitude');
  print('');
  print('🔬 STATISTICAL METHODOLOGY:');
  print('• Based on complex Wishart distribution test statistic');
  print('• Sequential testing identifies multiple change points');
  print('• Significance level controls false positive rate');
  print('• Median filtering reduces speckle noise effects');
  print('• Polarimetric information (VV+VH) enhances sensitivity');
  print('');
  print('📚 ALGORITHM REFERENCES:');
  print('• Nielsen, Conradsen & Skriver (2016) - Omnibus test methodology');
  print('• Conradsen et al. (2003) - Complex Wishart distribution theory');
  print('• Morton Canty - Google Earth Engine implementation');
  print('');
  print('✅ Analysis complete. Check Tasks tab for exports and click on map for details.');
  
} // End of main else block

// =========================================================================
// 10. Usage Instructions
// =========================================================================
print('');
print('🎛️ OMNIBUS ALGORITHM PARAMETERS:');
print('• significance: Statistical significance level (default: 0.0001)');
print('• relorbitnumber: Specific orbit number or "any" for all orbits');
print('• orbitpass: ASCENDING or DESCENDING orbit direction');
print('• median: Apply median filtering to reduce speckle noise');
print('• startDate/endDate: Analysis period (YYYY-MM-DD format)');
print('');
print('💡 OMNIBUS ALGORITHM TIPS:');
print('• Lower significance = fewer false alarms but may miss subtle changes');
print('• Higher significance = more sensitive but higher false alarm rate');
print('• Use specific orbit number to reduce computational load');
print('• Median filtering recommended for noisy data');
print('• Optimal time series length: 20-50 images');
print('• Requires dual-polarization data (VV+VH) for full effectiveness');
print('');
print('⚠️ COMPUTATIONAL CONSIDERATIONS:');
print('• Processing time scales quadratically with number of images');
print('• Large regions may require tiling for processing');
print('• Maximum recommended: 100 images, 50,000 km² region');
print('• Use shorter time periods or specific orbits to reduce load');