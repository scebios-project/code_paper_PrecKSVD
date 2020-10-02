clear all
close all

signal_numbers = [100:109 111:119 121:124 201:203 205 207:210 212:215 217 219:223 228 230:234];

disp('Reading samples ECG signal from MIT-BIH Arrhythmia Database')
wb = waitbar(0,'Reading samples ECG signal from MIT-BIH Arrhythmia Database');
for i=1:numel(signal_numbers)
    signal_name = num2str(signal_numbers(i));
    
    disp(['Reading signal ' signal_name])
    waitbar(i / numel(signal_numbers), wb, ['Reading signal ' signal_name]);
    
    [ecg_curr, Fs_curr, tm_curr]=rdsamp(['mitdb/' signal_name]);
    data{i} = ecg_curr;
    Fs{i}   = Fs_curr;
    tm{i}   = tm_curr;
    
    [ann_curr, type_curr, subtype_curr, chan_curr, num_curr]=rdann(['mitdb/' signal_name],'atr');
    ann{i}     = ann_curr;
    type{i}    = type_curr;
    subtype{i} = subtype_curr;
    chan{i}    = chan_curr;
    num{i}     = num_curr;
end
delete(wb);

save('data/MITDB.mat', 'data', 'Fs', 'tm', 'ann', 'type', 'subtype', 'chan', 'num');
%save('data/MITDB.mat', 'signal_numbers', '-append');