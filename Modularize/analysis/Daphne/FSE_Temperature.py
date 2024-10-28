import numpy as np

# Your data categorized by temperature
data_sets ={"15mK":"""
    25-10-24,18:50:50,1.477200e-02
    25-10-24,18:51:50,1.475700e-02
    25-10-24,18:52:50,1.476700e-02
    25-10-24,18:53:50,1.477700e-02
    25-10-24,18:54:50,1.477900e-02
    25-10-24,18:55:50,1.477100e-02
    25-10-24,18:56:50,1.476200e-02
    25-10-24,18:57:50,1.475900e-02
    25-10-24,18:58:50,1.476000e-02
    25-10-24,18:59:50,1.476700e-02
    25-10-24,19:00:50,1.475600e-02
    25-10-24,19:01:50,1.476100e-02
    25-10-24,19:02:50,1.474200e-02
    25-10-24,19:03:50,1.476900e-02
    25-10-24,19:04:50,1.473200e-02
    25-10-24,19:05:50,1.475300e-02
    25-10-24,19:06:50,1.476400e-02
    25-10-24,19:07:50,1.476500e-02
    25-10-24,19:08:50,1.474300e-02
    25-10-24,19:09:50,1.476700e-02
    25-10-24,19:10:50,1.475800e-02
    25-10-24,19:11:50,1.475300e-02
    25-10-24,19:12:50,1.474600e-02
    25-10-24,19:13:50,1.476000e-02
    """,
    "20mK": """
    26-10-24,18:25:50,2.051900e-02
    26-10-24,18:26:50,2.054800e-02
    26-10-24,18:27:50,2.060300e-02
    26-10-24,18:28:50,2.062600e-02
    26-10-24,18:29:50,2.066600e-02
    26-10-24,18:30:50,2.068200e-02
    26-10-24,18:31:50,2.069600e-02
    26-10-24,18:32:50,2.072400e-02
    26-10-24,18:33:50,2.076800e-02
    26-10-24,18:34:50,2.078600e-02
    26-10-24,18:35:50,2.078500e-02
    26-10-24,18:36:50,2.079800e-02
    26-10-24,18:37:50,2.083700e-02
    26-10-24,18:38:50,2.084400e-02
    26-10-24,18:39:50,2.085000e-02
    26-10-24,18:40:50,2.084100e-02
    """,
    "25mK": """
    26-10-24,21:35:50,2.495600e-02
    26-10-24,21:36:50,2.520300e-02
    26-10-24,21:37:50,2.522900e-02
    26-10-24,21:38:50,2.534600e-02
    26-10-24,21:39:50,2.535100e-02
    26-10-24,21:40:50,2.538000e-02
    26-10-24,21:41:50,2.535200e-02
    26-10-24,21:42:50,2.541000e-02
    26-10-24,21:43:50,2.547500e-02
    26-10-24,21:44:50,2.550900e-02
    26-10-24,21:45:50,2.544500e-02
    26-10-24,21:46:50,2.543900e-02
    26-10-24,21:47:50,2.547300e-02
    26-10-24,21:48:50,2.547500e-02
    26-10-24,21:49:50,2.550400e-02
    26-10-24,21:50:50,2.546600e-02
    26-10-24,21:51:50,2.546100e-02
    26-10-24,21:52:50,2.547000e-02
    26-10-24,21:53:50,2.554200e-02
    26-10-24,21:54:50,2.557300e-02
    26-10-24,21:55:50,2.550500e-02
    26-10-24,21:56:50,2.544200e-02
    26-10-24,21:57:50,2.544200e-02
    26-10-24,21:58:50,2.556700e-02
    26-10-24,21:59:50,2.551900e-02
    26-10-24,22:00:50,2.549000e-02
    26-10-24,22:01:50,2.554700e-02
    26-10-24,22:02:50,2.551300e-02
    26-10-24,22:03:50,2.553900e-02
    26-10-24,22:04:50,2.555100e-02
    26-10-24,22:05:50,2.549800e-02 
    """,
    "30mK": """
    26-10-24,23:44:50,3.047400e-02
    26-10-24,23:45:50,3.053100e-02
    26-10-24,23:46:50,3.048900e-02
    26-10-24,23:47:50,3.054500e-02
    26-10-24,23:48:50,3.047600e-02
    26-10-24,23:49:50,3.045300e-02
    26-10-24,23:50:50,3.045600e-02
    26-10-24,23:51:50,3.057800e-02
    26-10-24,23:52:50,3.053400e-02
    26-10-24,23:53:50,3.054300e-02
    26-10-24,23:54:50,3.058600e-02
    26-10-24,23:55:50,3.058800e-02
    26-10-24,23:56:50,3.056300e-02
    26-10-24,23:57:50,3.053700e-02
    26-10-24,23:58:50,3.053700e-02
    26-10-24,23:59:50,3.059700e-02
    27-10-24,00:00:50,3.061900e-02
    27-10-24,00:01:50,3.062000e-02
    27-10-24,00:02:50,3.059000e-02
    27-10-24,00:03:50,3.048700e-02
    """,
    "35mK": """
    27-10-24,12:05:50,3.548400e-02
    27-10-24,12:06:50,3.546900e-02
    27-10-24,12:07:50,3.556900e-02
    27-10-24,12:08:50,3.541000e-02
    27-10-24,12:09:50,3.550300e-02
    27-10-24,12:10:50,3.557100e-02
    27-10-24,12:11:50,3.556700e-02
    27-10-24,12:12:50,3.552500e-02
    27-10-24,12:13:50,3.556100e-02
    27-10-24,12:14:50,3.554200e-02
    27-10-24,12:15:50,3.557900e-02
    27-10-24,12:16:50,3.558500e-02
    27-10-24,12:17:50,3.562500e-02
    27-10-24,12:18:50,3.575800e-02
    27-10-24,12:19:50,3.562000e-02
    27-10-24,12:20:50,3.574800e-02
    27-10-24,12:21:50,3.580400e-02
    27-10-24,12:22:50,3.576400e-02
    27-10-24,12:23:50,3.578700e-02
    27-10-24,12:24:50,3.582300e-02
    """,
    "40mK": """
    27-10-24,14:20:50,4.042100e-02
    27-10-24,14:21:50,4.045800e-02
    27-10-24,14:22:50,4.041900e-02
    27-10-24,14:23:50,4.042600e-02
    27-10-24,14:24:50,4.035600e-02
    27-10-24,14:25:50,4.045400e-02
    27-10-24,14:26:50,4.053500e-02
    27-10-24,14:27:50,4.044000e-02
    27-10-24,14:28:50,4.041100e-02
    27-10-24,14:29:50,4.045300e-02
    27-10-24,14:30:50,4.042100e-02
    27-10-24,14:31:50,4.052700e-02
    27-10-24,14:32:50,4.034400e-02
    27-10-24,14:33:50,4.045900e-02
    27-10-24,14:34:50,4.033500e-02
    27-10-24,14:35:50,4.047000e-02
    27-10-24,14:36:50,4.041700e-02
    27-10-24,14:37:50,4.054600e-02
    27-10-24,14:38:50,4.051900e-02
    27-10-24,14:39:50,4.055600e-02
    27-10-24,14:40:50,4.037900e-02
    27-10-24,14:41:50,4.062900e-02
    27-10-24,14:42:50,4.049100e-02
    """,
    "45mK": """
    27-10-24,16:19:50,4.513100e-02
    27-10-24,16:20:50,4.508800e-02
    27-10-24,16:21:50,4.510100e-02
    27-10-24,16:22:50,4.507200e-02
    27-10-24,16:23:50,4.516000e-02
    27-10-24,16:24:50,4.515200e-02
    27-10-24,16:25:50,4.510500e-02
    27-10-24,16:26:50,4.510700e-02
    27-10-24,16:27:50,4.509700e-02
    27-10-24,16:28:50,4.513800e-02
    27-10-24,16:29:50,4.512300e-02
    27-10-24,16:30:50,4.513200e-02
    27-10-24,16:31:50,4.511700e-02
    27-10-24,16:32:50,4.509500e-02
    27-10-24,16:33:50,4.517800e-02
    27-10-24,16:34:50,4.512600e-02
    27-10-24,16:35:50,4.516400e-02
    27-10-24,16:36:50,4.511300e-02
    27-10-24,16:37:50,4.514600e-02
    27-10-24,16:38:50,4.521200e-02
    27-10-24,16:39:50,4.520300e-02
    """,
    "50mK": """
    27-10-24,18:38:50,4.991700e-02
    27-10-24,18:39:50,4.993800e-02
    27-10-24,18:40:50,4.985700e-02
    27-10-24,18:41:50,4.988800e-02
    27-10-24,18:42:50,4.991300e-02
    27-10-24,18:43:50,4.996600e-02
    27-10-24,18:44:50,4.991500e-02
    27-10-24,18:45:50,4.990200e-02
    27-10-24,18:46:50,4.992200e-02
    27-10-24,18:47:50,4.994900e-02
    27-10-24,18:48:50,4.993300e-02
    27-10-24,18:49:50,4.992300e-02
    27-10-24,18:50:50,4.989600e-02
    27-10-24,18:51:50,4.993400e-02
    27-10-24,18:52:50,4.991700e-02
    27-10-24,18:53:50,4.984900e-02
    27-10-24,18:54:50,4.986000e-02
    27-10-24,18:55:50,4.984800e-02
    27-10-24,18:56:50,4.988800e-02
    27-10-24,18:57:50,4.990800e-02
    """,
    "55mK": """
    27-10-24,21:22:50,5.535700e-02
    27-10-24,21:23:50,5.539700e-02
    27-10-24,21:24:50,5.531100e-02
    27-10-24,21:25:50,5.539800e-02
    27-10-24,21:26:50,5.534600e-02
    27-10-24,21:27:50,5.538400e-02
    27-10-24,21:28:50,5.538400e-02
    27-10-24,21:29:50,5.536500e-02
    27-10-24,21:30:50,5.535300e-02
    27-10-24,21:31:50,5.539900e-02
    27-10-24,21:32:50,5.538400e-02
    27-10-24,21:33:50,5.539600e-02
    27-10-24,21:34:50,5.543200e-02
    27-10-24,21:35:50,5.543500e-02
    27-10-24,21:36:50,5.545500e-02
    27-10-24,21:37:50,5.550900e-02
    27-10-24,21:38:50,5.549500e-02
    27-10-24,21:39:50,5.541200e-02
    27-10-24,21:40:50,5.544600e-02
    27-10-24,21:41:50,5.539100e-02
    27-10-24,21:42:50,5.548800e-02
    27-10-24,21:43:50,5.547100e-02 
    """,
    "60mK": """
    27-10-24,23:29:50,6.009200e-02
    27-10-24,23:30:50,6.009000e-02
    27-10-24,23:31:50,6.007900e-02
    27-10-24,23:32:50,6.014600e-02
    27-10-24,23:33:50,6.010000e-02
    27-10-24,23:34:50,6.012300e-02
    27-10-24,23:35:50,6.005000e-02
    27-10-24,23:36:50,6.018000e-02
    27-10-24,23:37:50,6.012800e-02
    27-10-24,23:38:50,6.022200e-02
    27-10-24,23:39:50,6.007200e-02
    27-10-24,23:40:50,6.012500e-02
    27-10-24,23:41:50,6.017400e-02
    27-10-24,23:42:50,6.018000e-02
    27-10-24,23:43:50,6.008500e-02
    27-10-24,23:44:50,6.012800e-02
    27-10-24,23:45:50,6.013900e-02
    27-10-24,23:46:50,6.017600e-02
    27-10-24,23:47:50,6.019400e-02
    27-10-24,23:48:50,6.020200e-02
    27-10-24,23:49:50,6.015200e-02
    """
} 

# Function to process each dataset
def process_dataset(data_str):
    lines = data_str.strip().split('\n')
    values = [float(line.split(',')[2]) * 1e3 for line in lines]
    values_array = np.array(values)
    mean = np.mean(values_array)
    std_dev = np.std(values_array)
    return mean, std_dev

# Process and print mean and standard deviations for all datasets
for temp, data_str in data_sets.items():
    mean, std_dev = process_dataset(data_str)
    print(f"{temp} Mean: {mean:.3f} +/- Std Dev: {std_dev:.4f}")