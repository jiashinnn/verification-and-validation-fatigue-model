"""
Questionnaire Data Processing Script
=====================================
Converts raw questionnaire data to normalized 0-1 range for fatigue simulation model.

Input:  Questionnaire_Data.csv (raw data)
Output: Questionnaire_Data_Processed.csv (normalized data)

Author: Fatigue Simulation Research Team
Date: 2025-12-26
"""

import pandas as pd
import numpy as np


# ============================================================================
# CONVERSION FUNCTIONS
# ============================================================================

def convert_heart_rate(hr_before, hr_after):
    """
    Convert heart rate measurements to 0-1 range based on gap.
    Uses discrete ranges matching backend PHP logic.

    Args:
        hr_before: Heart rate before training (bpm)
        hr_after: Heart rate after training (bpm)

    Returns:
        Hr value in [0, 1]
    """
    gap = float(hr_after) - float(hr_before)

    if gap <= 0:
        return 0.0
    if gap <= 10:
        return 0.2
    if gap <= 20:
        return 0.5
    if gap <= 30:
        return 0.75
    if gap <= 40:
        return 0.9
    return 1.0


def convert_sleep_hours(hours):
    """
    Convert sleep hours to 0-1 range.
    Optimal zone: 7-8 hours = 1.0
    Below 7h: faster decrease (0.15 per hour)
    Above 8h: slower decrease (0.05 per hour)

    Args:
        hours: Sleep hours (0-24)

    Returns:
        Sh value in [0, 1]
    """
    h = float(hours)

    if h <= 0.0:
        return 0.0

    # Perfect zone: 7-8 hours
    if 7.0 <= h <= 8.0:
        return 1.0

    # Below 7h: faster decrease (0.15 per hour)
    if h < 7.0:
        missing = 7.0 - h
        return max(0.0, min(1.0, 1.0 - (missing * 0.15)))

    # Above 8h: slower decrease (0.05 per hour)
    extra = h - 8.0
    return max(0.0, min(1.0, 1.0 - (extra * 0.05)))


def convert_sleep_timing(time_str):
    """
    Convert sleep timing to 0-1 range using piecewise linear interpolation.
    Based on anchor points: 22:30→1.0, 01:00→0.5, 02:00→0.4, ..., >05:00→0.0

    Args:
        time_str: Time in "HH:MM:SS" format (e.g., "22:30:00", "02:00:00")

    Returns:
        St value in [0, 1]
    """
    # Parse time string
    parts = time_str.strip().split(':')
    hh = int(parts[0])
    mm = int(parts[1])

    # Convert to decimal hours since midnight
    t = hh + (mm / 60.0)

    # Treat after-midnight (00:00-11:59) as next day (24:00-35:59)
    if t < 12.0:
        t += 24.0

    # Anchor points (in hours since midnight, adjusted for next-day)
    ideal = 22.5    # 22:30 → 1.0
    a1 = 25.0       # 01:00 → 0.5
    a2 = 26.0       # 02:00 → 0.4
    a3 = 27.0       # 03:00 → 0.3
    a4 = 28.0       # 04:00 → 0.2
    a5 = 29.0       # 05:00 → 0.1
    cutoff = 29.0   # After 05:00 → 0.0

    # Case 1: 22:30 or earlier → 1.0 (flat ceiling)
    if t <= ideal:
        return 1.0

    # Case 2: After 05:00 → 0.0 (too late)
    if t > cutoff:
        return 0.0

    # Case 3-7: Piecewise linear interpolation between anchor points
    if t <= a1:
        # Between 22:30 and 01:00: linear 1.0 → 0.5
        return 1.0 - ((t - ideal) / (a1 - ideal)) * (1.0 - 0.5)

    if t <= a2:
        # Between 01:00 and 02:00: linear 0.5 → 0.4
        return 0.5 - ((t - a1) / (a2 - a1)) * (0.5 - 0.4)

    if t <= a3:
        # Between 02:00 and 03:00: linear 0.4 → 0.3
        return 0.4 - ((t - a2) / (a3 - a2)) * (0.4 - 0.3)

    if t <= a4:
        # Between 03:00 and 04:00: linear 0.3 → 0.2
        return 0.3 - ((t - a3) / (a4 - a3)) * (0.3 - 0.2)

    if t <= a5:
        # Between 04:00 and 05:00: linear 0.2 → 0.1
        return 0.2 - ((t - a4) / (a5 - a4)) * (0.2 - 0.1)

    # Fallback (should not reach here due to cutoff check)
    return 0.0


def convert_likert(value):
    """
    Convert Likert scale (1-5) to 0-1 range.

    Args:
        value: Likert scale value (1-5)

    Returns:
        Normalized value in [0, 1]
    """
    likert_mapping = {
        1: 0.0,
        2: 0.4,
        3: 0.6,
        4: 0.8,
        5: 1.0
    }
    return likert_mapping.get(value, 0.0)


def calculate_fatigue(row):
    """
    Calculate fatigue score from 10 FAS (Fatigue Assessment Scale) questions.
    Items 4 and 10 are reverse scored.

    Args:
        row: DataFrame row containing FAS question columns

    Returns:
        Fatigue value in [0, 1]
    """
    # Column names for FAS questions (in order)
    fas_columns = [
        "I am bothered by fatigue",                              # Q1
        "I get tired very quickly",                              # Q2
        "I don't do much during the day",                        # Q3
        "I have enough energy for everyday life",                # Q4 - REVERSE
        "Physically, I feel exhausted",                          # Q5
        "I have problems starting things",                       # Q6
        "I have problems thinking clearly",                      # Q7
        "I feel no desire to do anything",                       # Q8
        "Mentally, I feel exhausted",                            # Q9
        "When I am doing something, I can concentrate quite well" # Q10 - REVERSE
    ]

    # Extract values
    scores = []
    for i, col in enumerate(fas_columns):
        value = row[col]

        # Reverse scoring for items 4 and 10 (indices 3 and 9)
        if i == 3 or i == 9:
            # Reverse: 5→1, 4→2, 3→3, 2→4, 1→5
            reversed_value = 6 - value
            scores.append(reversed_value)
        else:
            scores.append(value)

    # Sum all scores (range: 10-50)
    raw_score = sum(scores)

    # Normalize to 0-1 range
    # Formula: (raw_score - 10) / 40
    fatigue_normalized = (raw_score - 10) / 40.0

    return fatigue_normalized


# ============================================================================
# MAIN PROCESSING FUNCTION
# ============================================================================

def process_questionnaire_data(input_file, output_file):
    """
    Process raw questionnaire data and save normalized results.

    Args:
        input_file: Path to input CSV file
        output_file: Path to output CSV file
    """
    print("=" * 70)
    print("QUESTIONNAIRE DATA PROCESSING")
    print("=" * 70)

    # Read raw data
    print(f"\nReading data from: {input_file}")
    df = pd.read_csv(input_file)
    print(f"[OK] Loaded {len(df)} rows")

    # Create output dataframe with metadata columns
    print("\nProcessing data...")
    output_df = pd.DataFrame()

    # Copy metadata columns (first 5 columns)
    metadata_cols = ['Name', 'Service_Category', 'Level_of_Training', 'Date_of_Training', 'Gender']
    for col in metadata_cols:
        output_df[col] = df[col]

    # Convert Heart Rate (Hr)
    print("  [OK] Converting Heart Rate (Hr)...")
    output_df['Hr'] = df.apply(
        lambda row: convert_heart_rate(row['Heart_Rate_Before'], row['Heart_Rate_After']),
        axis=1
    ).round(2)

    # Convert Sleep Hours (Sh)
    print("  [OK] Converting Sleep Hours (Sh)...")
    output_df['Sh'] = df['Sleep_Hours'].apply(convert_sleep_hours).round(2)

    # Convert Sleep Timing (St)
    print("  [OK] Converting Sleep Timing (St)...")
    output_df['St'] = df['Sleep_Timing'].apply(convert_sleep_timing).round(2)

    # Convert Likert scale inputs
    print("  [OK] Converting Likert scale inputs (Tp, Fi, Fl, Te, Ti, Tc)...")
    output_df['Tp'] = df['Training_Pressure'].apply(convert_likert).round(2)
    output_df['Fi'] = df['Fluid_Intake'].apply(convert_likert).round(2)
    output_df['Fl'] = df['Fitness_Level'].apply(convert_likert).round(2)
    output_df['Te'] = df['Thermal_Environment'].apply(convert_likert).round(2)
    output_df['Ti'] = df['Training_Intensity'].apply(convert_likert).round(2)
    output_df['Tc'] = df['Training_Complexity'].apply(convert_likert).round(2)

    # Calculate Fatigue from FAS questions
    print("  [OK] Calculating Fatigue from FAS questionnaire...")
    output_df['Fatigue'] = df.apply(calculate_fatigue, axis=1).round(2)

    # Save processed data
    print(f"\nSaving processed data to: {output_file}")
    output_df.to_csv(output_file, index=False)
    print(f"[OK] Saved {len(output_df)} rows with {len(output_df.columns)} columns")

    # Print summary statistics
    print("\n" + "=" * 70)
    print("SUMMARY STATISTICS (Normalized Values)")
    print("=" * 70)

    input_cols = ['Hr', 'Sh', 'St', 'Tp', 'Fi', 'Fl', 'Te', 'Ti', 'Tc', 'Fatigue']

    print(f"\n{'Variable':<12} {'Min':<8} {'Max':<8} {'Mean':<8} {'Std':<8}")
    print("-" * 50)

    for col in input_cols:
        col_data = output_df[col]
        print(f"{col:<12} {col_data.min():<8.2f} {col_data.max():<8.2f} "
              f"{col_data.mean():<8.2f} {col_data.std():<8.2f}")

    print("\n" + "=" * 70)
    print("PROCESSING COMPLETE!")
    print("=" * 70)
    print(f"\nOutput file: {output_file}")
    print("\nTo add new data:")
    print("  1. Add new rows to 'Questionnaire_Data.csv'")
    print("  2. Run this script again: python process_questionnaire_data.py")
    print("  3. The processed file will be regenerated with all data")
    print("=" * 70)


# ============================================================================
# SCRIPT ENTRY POINT
# ============================================================================

if __name__ == "__main__":
    # File paths
    INPUT_FILE = "Questionnaire_Data.csv"
    OUTPUT_FILE = "Questionnaire_Data_Processed.csv"

    # Process data
    try:
        process_questionnaire_data(INPUT_FILE, OUTPUT_FILE)
    except FileNotFoundError:
        print(f"\n[ERROR] File '{INPUT_FILE}' not found!")
        print("Please ensure the file is in the same directory as this script.")
    except Exception as e:
        print(f"\n[ERROR] {str(e)}")
        import traceback
        traceback.print_exc()
