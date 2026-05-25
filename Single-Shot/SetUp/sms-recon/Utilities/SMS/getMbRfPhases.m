function [ rfPhases] = getMbRfPhases(mbFactor)
   %getMbRFPhases Convenience function to return phases for constructing SMS
   %pulses using optimized phase schedules to minimize peak power.
   % See Wong, E. "Optimized phase schedules for minimizing peak RF power in
   % simultaneous multi-slice RF excitation pulses." ISMRM 2012
   %
   %
   % Inputs
   %    mbFactor - Integer (1,9) that describes number of bands
   %
   % Outputs
   %    rfPhases - Row vector of real phases applied to each slice in the RF
   %               excitation. Length is [1,mbFactor];
   
   switch mbFactor
      case 1
         rfPhases = 0;
      case 2
         rfPhases = [0, 0];
      case 3
         rfPhases = [0, 0.730, 4.602];
      case 4
         rfPhases = [0, 3.875, 5.940, 6.197];
      case 5
         rfPhases = [0, 3.778, 5.335, 0.872, 0.471];
      case 6
         rfPhases = [0, 2.005, 1.674, 5.012, 5.736, 4.123];
      case 7
         rfPhases = [0, 3.002, 5.998, 5.909, 2.624, 2.528, 2.440];
      case 8
         rfPhases = [0, 1.036, 3.414, 3.778, 3.215, 1.756, 4.555, 2.467];
      case 9
         rfPhases = [0, 1.250, 1.783, 3.558, 0.739, 3.319, 1.296, 0.521, 0.5332];
      case 10
         rfPhases = [0, 4.418, 2.360, 0.677, 2.253, 3.472, 3.040, 3.974, 1.192, 2.510];
      case 11
         rfPhases = [0, 5.041, 4.285, 3.001, 5.765, 4.295, 0.056, 4.213, 6.040, 1.078, 2.759];
      case 12
         rfPhases = [0, 2.755, 5.491, 4.447, 0.231, 2.499, 3.539, 2.931, 2.759, 5.376, 4.554, 3.479];
      case 13
         rfPhases = [0, 0.603, 0.009, 4.179, 4.361, 4.837, 0.816, 5.995, 4.150, 0.417, 1.520, 4.517, 1.729];
      case 14
         rfPhases = [0, 3.997, 0.830, 5.712, 3.838, 0.084, 1.685, 5.328, 0.237, 0.506, 1.356, 4.025, 4.483, 4.084];
      case 15
         rfPhases = [0, 4.126, 2.266, 0.957, 4.603, 0.815, 3.475, 0.997, 1.449, 1.192, 0.148, 0.939, 2.531, 3.612, 4.801];
      case 16
         rfPhases = [0, 4.359, 3.510, 4.410, 1.750, 3.357, 2.061, 5.948, 3.000, 2.822, 0.627, 2.768, 3.875, 4.173, 4.224, 5.941];
      otherwise
         error('Bad Multiband Factor: Unsupported number of bands, only 1 - 16 bands are supported! Seriously, how many bands can you cram in there?');
   end
end

