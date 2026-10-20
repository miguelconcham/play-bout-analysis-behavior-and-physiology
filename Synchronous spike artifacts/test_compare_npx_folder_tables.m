function test_compare_npx_folder_tables()
%TEST_COMPARE_NPX_FOLDER_TABLES  Folder-composition checks on synthetic sheets.

root = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM';

info = {
    'Probe', 'Animal code', 'Recording number', 'Date', 'NPX recording folder'
    'ProbeA', 'B2S2', 1, 20231110, 'PlayRec'
    1, 'B1S3', 1, '2023-10-08', 'catgt'
    'ProbeB', 'B4D4', 2, 826, 'PlayRec'
    };
folders = {
    'Session name', 'Continuous file'
    'B2S2 1110 Single2', [root '\20231110\PlayRec\Neuropix-PXI-100.ProbeA-AP\recording1\continuous.dat']
    'B1S3 1008 Single', [root '/2023-10-08/catgt/catgt_g0/catgt_g0_imec0/recording1/continuous.dat']
    'B4D4 0826 Dual', [root '\20230826\PlayRec\recording2\Neuropix-PXI-100.ProbeB-AP\continuous.dat']
    };
report = compare_npx_folder_tables(info, folders);
assert(report.ok, 'Matching sheets should agree.');
assert(report.n_experiment == 3, 'Empty trailing fields should not add rows.');
assert(report.n_bad == 0, 'A matching sheet should have no problems.');

% A second file for the same session is still a match, not an extra error.
folders_two = [folders; {
    'B2S2 1110 Single2', [root '\20231110\PlayRec\recording1\ProbeA\continuous.dat']
    }];
report = compare_npx_folder_tables(info, folders_two);
assert(report.ok, 'Two valid files for one session should agree.');

% Column order when the headers are not names we know.
ordered = {
    'A', 'B', 'C', 'D', 'E'
    'B3D2', '20231130', 'RecFolder', 'ProbeB', 2
    };
ordered_folders = {
    'Session', 'Path'
    'B3D2 1130 Dual', [root '\20231130\RecFolder\recording2\ProbeB\continuous.dat']
    };
report = compare_npx_folder_tables(ordered, ordered_folders);
assert(report.ok, 'Columns 1-5 should be used when headers are unrecognized.');
assert(~isempty(report.notes), 'Column-order fallback should be reported.');

% Excel serial for 2023-11-10, and an MMDD stored without the leading zero.
serial = datenum(2023, 11, 10) - datenum(1899, 12, 30);
serial_info = {
    'Animal', 'Date', 'RecordingFolder', 'Probe', 'RecordingNumber'
    'B4S2', serial, 'PlayRec', 1, 1
    'B4S2', 825, 'OtherRec', 'imec1', 'recording4'
    };
serial_folders = {
    'Session', 'Path'
    'B4S2 1110 Single', [root '\20231110\PlayRec\recording1\imec0\continuous.dat']
    'B4S2 0825 Single', [root '\20230825\OtherRec\imec1\recording4\continuous.dat']
    };
report = compare_npx_folder_tables(serial_info, serial_folders);
assert(report.ok, 'Excel serial dates and unpadded MMDD values should agree.');

bad_number = folders;
bad_number{2, 2} = [root '\20231110\PlayRec\Neuropix-PXI-100.ProbeA-AP\recording12\continuous.dat'];
report = compare_npx_folder_tables(info(1:2, :), bad_number(1:2, :));
assert(~report.ok, 'recording12 must not satisfy recording number 1.');
assert_problem(report, 'bad path', 'recording-number');

bad_probe = folders;
bad_probe{3, 2} = [root '/2023-10-08/catgt/catgt_g0/catgt_g0_imec1/recording1/continuous.dat'];
report = compare_npx_folder_tables(info([1 3], :), bad_probe([1 3], :));
assert(~report.ok, 'imec1 must not satisfy probe number 1.');
assert_problem(report, 'bad path', 'probe');

swapped = {
    'Session', 'Path'
    'B2S2 1110 Single2', [root '\PlayRec\20231110\recording1\ProbeA\continuous.dat']
    };
report = compare_npx_folder_tables(info(1:2, :), swapped);
assert_problem(report, 'bad path', 'NPX recording folder');

animal_after = {
    'Session', 'Path'
    'B2S2 1110 Single2', [root '\20231110\B2S2\PlayRec\recording1\ProbeA\continuous.dat']
    };
report = compare_npx_folder_tables(info(1:2, :), animal_after);
assert_problem(report, 'bad path', 'after the date');

wrong_day = {
    'Session', 'Path'
    'B2S2 1111 Single2', [root '\20231111\PlayRec\recording1\ProbeA\continuous.dat']
    };
report = compare_npx_folder_tables(info(1:2, :), wrong_day);
assert_problem(report, 'no session', 'B2S2 1110');
assert_problem(report, 'unlisted', 'not in the experiment info');

wrong_label = {
    'Session', 'Path'
    'B2S2 1110 Dual', [root '\20231110\PlayRec\recording1\ProbeA\continuous.dat']
    };
report = compare_npx_folder_tables(info(1:2, :), wrong_label);
assert_problem(report, 'bad session name', 'Dual');

unparsed = {
    'Session', 'Path'
    'not-a-session', [root '\20231110\PlayRec\recording1\ProbeA\continuous.dat']
    };
report = compare_npx_folder_tables(info(1:2, :), unparsed);
assert_problem(report, 'bad session name', 'not-a-session');
assert_problem(report, 'no session', 'B2S2 1110');

area_limits = {
    'AnimalName', 'depth_start', 'depth_end', 'ProbeNum'
    'Batch2Single2', 1000, 2000, 1
    };
threw = false;
caught = '';
try
    compare_npx_folder_tables(area_limits, folders);
catch ME
    threw = true;
    caught = ME.identifier;
end
assert(threw, 'The depth workbook should be rejected.');
assert(strcmp(caught, 'compare_npx_folder_tables:areaLimits'), ...
    ['Expected areaLimits, got ', caught, '.']);

empty_info = {'Animal', 'Date', 'RecordingFolder', 'Probe', 'RecordingNumber'};
empty_folders = {'Session', 'Path'};
report = compare_npx_folder_tables(empty_info, empty_folders);
assert(~report.ok, 'An empty experiment sheet should not count as agreement.');
assert(report.n_bad == 0, 'An empty sheet has no mismatch rows.');

fprintf('compare_npx_folder_tables tests passed.\n');
end

function assert_problem(report, status, fragment)
hit = false;
for k = 1:numel(report.rows)
    row = report.rows(k);
    if strcmp(row.status, status) && ~isempty(strfind(row.detail, fragment))
        hit = true;
    end
end
assert(hit, ['Missing ' status ' containing "' fragment '".']);
end
