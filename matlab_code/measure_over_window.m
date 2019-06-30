function measurement_over_window = measure_over_window(arclength_and_measurement,smoothed_line_delta, marker, window_in_mm)
%{
INPUTS:
    arclength_and_measurement: [nx2] array, where col1 is arclength (eg [0;
    0.1; 0.2; 0.3; ... n]), and col2 is the measurement at that given
    arclength. 

    smoothed_line_delta: [mx4] array, where col1 is arclength and corresponds
    exactly to the arclength_and_measurement variable. 

    marker: a point on smoothed_line_delta from which we wish to look at a
    specific measurement.

    window_in_mm: a window to take all measurements, starting a specific
    marker point and then ending after the windowed distance. 


%}
smoothed_line = smoothed_line_delta(:,2:4);
window_in_10mm = window_in_mm * 10; % increments are in 10 units per mm.

markerIndexOnLine = find_index_of_row_in_array(marker,smoothed_line);
markerOnLine = round(smoothed_line_delta(markerIndexOnLine,1),1);

curve_at_marker = arclength_and_measurement(round(arclength_and_measurement(:,1),1)==markerOnLine,:);
curve_at_marker_index = find_index_of_row_in_array(curve_at_marker,arclength_and_measurement);

curves_at_marker_over_window = arclength_and_measurement(curve_at_marker_index:curve_at_marker_index+window_in_10mm,:);
measurement_over_window = {mean(curves_at_marker_over_window(:,2)),curves_at_marker_over_window};

end