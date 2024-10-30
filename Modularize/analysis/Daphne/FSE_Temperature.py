import os
import re
import numpy as np
from datetime import datetime

# Find the earliest and latest times from filenames in a folder
def find_time_range_from_filenames(folder_path, log_date):
    times = []
    pattern = re.compile(r"H(\d{1,2})M(\d{2})S(\d{2})")  # Match single or double digit hour format

    for filename in os.listdir(folder_path):
        if filename.endswith(".nc"):  # Only process .nc files
            match = pattern.search(filename)
            if match:
                # Extract time and add a leading zero if hour is a single digit
                hour, minute, second = match.groups()
                hour = hour.zfill(2)  # Ensure hour has two digits
                time_str = f"{log_date},{hour}:{minute}:{second}"
                time_obj = datetime.strptime(time_str, '%y-%m-%d,%H:%M:%S')
                times.append(time_obj)
    
    # Find the earliest and latest times
    if times:
        start_time = min(times)
        end_time = max(times)
        return start_time.strftime('%y-%m-%d,%H:%M:%S'), end_time.strftime('%y-%m-%d,%H:%M:%S')
    else:
        return None, None

# Read data from a .log file and calculate the mean and standard deviation
def process_file(file_path, start_time_str, end_time_str):
    start_time = datetime.strptime(start_time_str, '%y-%m-%d,%H:%M:%S')
    end_time = datetime.strptime(end_time_str, '%y-%m-%d,%H:%M:%S')
    values = []

    with open(file_path, 'r') as file:
        for line in file:
            date_str, time_str, value_str = line.strip().split(',')
            timestamp = datetime.strptime(f"{date_str},{time_str}", '%d-%m-%y,%H:%M:%S')
            value = float(value_str) * 1e3
            
            if start_time <= timestamp <= end_time:
                values.append(value)

    values_array = np.array(values)
    if values_array.size > 0:
        mean = np.mean(values_array)
        std_dev = np.std(values_array)
        return mean, std_dev
    else:
        print("No data found in the specified time range")
        return None, None

# Recursively find folders that match target names within the main folder
def find_target_folders(main_folder, target_folders):
    matched_folders = []
    for root, dirs, files in os.walk(main_folder):
        for dir_name in dirs:
            if dir_name in target_folders:
                matched_folders.append(os.path.join(root, dir_name))
    return matched_folders

# Main function to process each target folder's .nc files
def main(main_folder, log_file_path, target_folders_list):
    # Extract date from .log filename
    log_filename = os.path.basename(log_file_path)
    log_date_match = re.search(r"(\d{2}-\d{2}-\d{2})", log_filename)
    if log_date_match:
        log_date = log_date_match.group(1)  # Date format is YY-MM-DD
        print(f"Extracted date: {log_date}")

        # Specify the target folder names to search for
        target_folders = target_folders_list  # List of folder names to analyze
        matched_folders = find_target_folders(main_folder, target_folders)

        # Process each matched folder
        for folder in matched_folders:
            print(f"\nProcessing folder: {folder}")
            start_time_str, end_time_str = find_time_range_from_filenames(folder, log_date)
            if start_time_str and end_time_str:
                print(f"Time range: {start_time_str} to {end_time_str}")
                
                # Process .log file based on the time range
                mean, std_dev = process_file(log_file_path, start_time_str, end_time_str)
                if mean is not None:
                    print(f"Mean: {mean:.3f} +/- Std Dev: {std_dev:.4f}")
            else:
                print("No matching time found")
    else:
        print("Unable to extract date from filename")

# Define paths and target folders
main_folder = r"C:\Users\User\SynologyDrive\SynologyDrive\09 Data\Fridge Data\Qubit\20241024_DRKe_5XQv4#5_second_coating_and_effT\Meas_raw\Q3_CopyFoldersForMainAnalysis\QDbackupIs1025" # Main folder path
log_file_path = r"C:\Users\User\SynologyDrive\SynologyDrive\09 Data\Fridge Data\Qubit\20241024_DRKe_5XQv4#5_second_coating_and_effT\Meas_raw\Q3_CopyFoldersForMainAnalysis\QDbackupIs1025\CH9_FSE_24-10-25.log" # .log file name
target_folders_list = ['SS', 'T1', 'T2_Ramsey', 'T2_SpinEcho']

# Execute program
main(main_folder, log_file_path, target_folders_list)
