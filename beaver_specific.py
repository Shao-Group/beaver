#!/usr/bin/env python3

import argparse
import os
import csv
import math
import pandas as pd
import joblib
import numpy as np
from pathlib import Path

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Score transcripts using pre-trained Beaver-Specific model')
    parser.add_argument('-n', '--sample_size', type=int, required=True,
                      help='Number of cells in the dataset')
    parser.add_argument('-i', '--input', required=True, 
                      help='Input specific feature directory')
    parser.add_argument('-g', '--general_features', required=True,
                      help='CSV containing Beaver-General processed scores')
    parser.add_argument('-m', '--model', required=True,
                      help='Path to the pre-trained specific model file for scoring')
    parser.add_argument('-o', '--output', required=True,
                      help='Output directory')
    parser.add_argument('-p', '--probability', type=float, default=0.6,
                      help='Minimum specific probability score threshold (range: 0 to 1)')
    parser.add_argument('-pg', '--general_probability', type=float, default=0.2,
                      help='Minimum general probability score threshold for filtering (range: 0 to 1)')
    parser.add_argument('-s', '--save', action='store_true',
                      help='Save processed feature files and estimated scores')
    
    args = parser.parse_args()
    return args

def calc_stats(values):
    """Calculate statistics for a list of values."""
    if not values:
        return [0, 0, 0, 0, 0]  # Default values if the list is empty
    return [round(x, 3) for x in [np.min(values), np.max(values), np.median(values), np.mean(values), np.std(values)]]

def count_non_zero(values):
    """Count non-zero values in a list."""
    return sum(1 for v in values if v != 0)

def continuous_non_zero(values):
    """Calculate max and min continuous non-zero values."""
    if not values:
        return 0, 0
    
    counts = []
    current_count = 0
    
    for v in values:
        if v != 0:
            current_count += 1
        else:
            if current_count > 0:
                counts.append(current_count)
                current_count = 0
    
    if current_count > 0:
        counts.append(current_count)
    
    if not counts:
        return 0, 0
    
    return max(counts), min(counts)

def get_specific_feature_vector(v, sample_size):
    """Extract feature vector from raw data."""
    features = []
    
    # Basic information
    transcript_id, chrom = v[0], v[1]
    features.extend([transcript_id, chrom])
    
    # Input and output counts
    input_genes, input_transcripts, output_transcripts = int(v[2]), int(v[3]), int(v[4])
    features.extend([input_genes, input_transcripts, output_transcripts])
    
    # Transcript support information
    supporting_junctions = int(v[5])
    ratio_supporting_junctions = round(float(v[6]), 3)
    features.extend([supporting_junctions, ratio_supporting_junctions])
    # Number of junctions
    num_junctions = int(v[7])
    #features.append(num_junctions)
    
    # Process cell_comp_junc_cov
    start_idx = 8
    end_idx = start_idx + num_junctions
    cell_comp_junc_cov = [float(x) for x in v[start_idx:end_idx]]
    features.extend(calc_stats(cell_comp_junc_cov))
    
    # New features for cell_comp_junc_cov
    features.append(count_non_zero(cell_comp_junc_cov))
    max_continuous, min_continuous = continuous_non_zero(cell_comp_junc_cov)
    features.extend([max_continuous, min_continuous])
    
    # Process cell_junc_cov
    start_idx = end_idx
    end_idx = start_idx + num_junctions
    cell_junc_cov = [float(x) for x in v[start_idx:end_idx]]
    features.extend(calc_stats(cell_junc_cov))
    
    # New features for cell_junc_cov
    features.append(count_non_zero(cell_junc_cov))
    max_continuous, min_continuous = continuous_non_zero(cell_junc_cov)
    features.extend([max_continuous, min_continuous])
    
    return features

def process_cell_features(input_dir, sample_size, general_features, pg):
    """Process features for all cells and combine with general results."""
    all_data = []
    
    # Load general results and features
    general_df = pd.read_csv(general_features)
    general_df = general_df[general_df['general_prob'] >= pg]
    print(f"Number of transcripts passing general probability threshold ({pg}): {len(general_df)}")

    print(f"Processing features for {sample_size} cells...")
    for cell_id in range(1, sample_size + 1):
        feature_file = os.path.join(input_dir, f"{cell_id}_feature.csv")
        
        if not os.path.exists(feature_file):
            print(f"Warning: Feature file not found for cell {cell_id}")
            continue
            
        # Process feature file
        with open(feature_file, 'r') as f:
            reader = csv.reader(f)
            for row in reader:
                features = get_specific_feature_vector(row, sample_size)
                if features is None:
                    continue
                
                data_row = [cell_id] + features
                all_data.append(data_row)
    
    # Create DataFrame with cell-specific features
    feature_names = [
        'cell_id', 'transcript_id', 'chr',
        'input_genes', 'input_transcripts', 'output_transcripts',
        'cell_supporting_junctions', 'cell_ratio_supporting_junctions',
        'cell_min_comp_cov', 'cell_max_comp_cov', 'cell_median_comp_cov', 'cell_mean_comp_cov', 'cell_std_comp_cov',
        'cell_nonzero_comp_cov', 'cell_max_streak_comp_cov', 'cell_min_streak_comp_cov',
        'cell_min_junc_cov', 'cell_max_junc_cov', 'cell_median_junc_cov', 'cell_mean_junc_cov', 'cell_std_junc_cov',
        'cell_nonzero_junc_cov', 'cell_max_streak_junc_cov', 'cell_min_streak_junc_cov'
    ]
    
    df = pd.DataFrame(all_data, columns=feature_names)

    # Merge with general features
    general_feature_cols = [
        'transcript_id','sample_size',
        'bottleneck_coverage', 'highest_coverage',
        'extendable_score', 'num_junctions',
        'min_comp_cov', 'max_comp_cov', 'median_comp_cov', 'mean_comp_cov', 'std_comp_cov',
        'min_junc_cov', 'max_junc_cov', 'median_junc_cov', 'mean_junc_cov', 'std_junc_cov',
        'cells_support', 'ratio_cells_support',
        'min_cell_junc', 'max_cell_junc', 'median_cell_junc', 'mean_cell_junc', 'std_cell_junc',
        'ratio_cells_fl_over_sample', 'ratio_cells_fl_over_cells_support',
        'num_fragments',
        'min_frag_cov', 'max_frag_cov', 'median_frag_cov', 'mean_frag_cov', 'std_frag_cov',
        'general_prob'
    ]

    df = pd.merge(df, general_df[general_feature_cols], on='transcript_id', how='right')
    return df

def score_specific_transcripts(df, model, prob_threshold):
    """Score transcripts using the specific model."""

    feature_columns = [
        'sample_size',
        'input_genes', 'input_transcripts', 'output_transcripts',
        'cell_supporting_junctions', 'cell_ratio_supporting_junctions',
        'cell_min_comp_cov', 'cell_max_comp_cov', 'cell_median_comp_cov', 'cell_mean_comp_cov', 'cell_std_comp_cov',
        'cell_nonzero_comp_cov', 'cell_max_streak_comp_cov', 'cell_min_streak_comp_cov',
        'cell_min_junc_cov', 'cell_max_junc_cov', 'cell_median_junc_cov', 'cell_mean_junc_cov', 'cell_std_junc_cov',
        'cell_nonzero_junc_cov', 'cell_max_streak_junc_cov', 'cell_min_streak_junc_cov',
        'bottleneck_coverage', 'highest_coverage',
        'extendable_score', 'num_junctions',
        'min_comp_cov', 'max_comp_cov', 'median_comp_cov', 'mean_comp_cov', 'std_comp_cov',
        'min_junc_cov', 'max_junc_cov', 'median_junc_cov', 'mean_junc_cov', 'std_junc_cov',
        'cells_support', 'ratio_cells_support',
        'min_cell_junc', 'max_cell_junc', 'median_cell_junc', 'mean_cell_junc', 'std_cell_junc',
        'ratio_cells_fl_over_sample', 'ratio_cells_fl_over_cells_support',
        'num_fragments',
        'min_frag_cov', 'max_frag_cov', 'median_frag_cov', 'mean_frag_cov', 'std_frag_cov'
    ]
    
    # Get probability scores
    probabilities = model.predict_proba(df[feature_columns])[:, 1]
    
    # Filter results based on probability threshold
    mask = probabilities >= prob_threshold
    filtered_df = df[mask].copy()
    filtered_probabilities = probabilities[mask]
    
    # Create results DataFrame
    results = pd.DataFrame({
        'cell_id': filtered_df['cell_id'],
        'transcript_id': filtered_df['transcript_id'],
        'chr': filtered_df['chr'],
        'specific_probability': filtered_probabilities
    })
    
    return results

def main():
    args = parse_arguments()
    os.makedirs(args.output, exist_ok=True)
    
    # Load pre-trained model
    print(f"Loading pre-trained specific model from {args.model}")
    model = joblib.load(args.model)
    
    # Process cell features and combine with general results
    print("Processing cell-specific features...")
    df = process_cell_features(args.input, args.sample_size, args.general_features, args.general_probability)
    print(f"Processed features for {len(df)} transcript-cell combinations")
    
    # Save processed features
    if args.save:
        features_file = os.path.join(args.output, "specific_features.csv")
        df.to_csv(features_file, index=False)
        print(f"Saved processed specific features to {features_file}")
    
    # Score transcripts
    print("Scoring transcripts with specific model...")
    results = score_specific_transcripts(df, model, args.probability)
    
    # Save results
    output_file = os.path.join(args.output, "beaver_specific_scores.csv")
    results.to_csv(output_file, index=False)
    print(f"Saved scoring results to {output_file}")
    
    # Print summary statistics
    print("\nScoring Summary:")
    print(f"Total transcript-cell combinations: {len(df)}")
    print(f"Combinations passing threshold: {len(results)}")
    print(f"Combinations filtered out: {len(df) - len(results)}")
    print(f"Using general probability threshold: {args.general_probability}")
    print(f"Using specific probability threshold: {args.probability}")
    print(f"Unique transcripts in final results: {results['transcript_id'].nunique()}")
    print(f"Unique cells in final results: {results['cell_id'].nunique()}")

if __name__ == "__main__":
    main()