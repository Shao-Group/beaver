#!/usr/bin/env python3

import argparse
import os
import csv
import math
import pandas as pd
import joblib
import numpy as np

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Score transcripts using pre-trained Beaver-General model')
    parser.add_argument('-i', '--input', required=True, 
                      help='Input general feature CSV file')
    parser.add_argument('-m', '--model', required=True,
                      help='Path to the pre-trained model file for scoring')
    parser.add_argument('-n', '--sample_size', type=int, required=True,
                      help='Number of cells in the dataset')
    parser.add_argument('-o', '--output', required=True,
                      help='Output directory')
    parser.add_argument('-p', '--probability', type=float, default=0.2,
                      help='Minimum probability score threshold (range: 0 to 1)')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        raise ValueError(f"Input file does not exist: {args.input}")
    if not os.path.exists(args.model):
        raise ValueError(f"Model file does not exist: {args.model}")
    if args.probability < 0 or args.probability > 1:
        raise ValueError(f"Probability threshold must be between 0 and 1")
    if args.sample_size <= 0:
        raise ValueError(f"Sample size must be positive")
    
    return args

def calc_stats(values):
    """Calculate statistics for a list of values."""
    if not values:
        return [0, 0, 0, 0, 0]  # Default values if the list is empty
    return [round(x, 3) for x in [np.min(values), np.max(values), np.median(values), np.mean(values), np.std(values)]]

def get_meta_feature_vector(v, sample_size):
    """Extract feature vector from raw data."""
    try:
        features = []
        # Basic information
        transcript_id, chrom = v[0], v[1]
        features.extend([transcript_id, chrom, sample_size])

        # Transcript support information
        bottleneck_coverage = round(float(v[2]), 3)
        highest_coverage = round(float(v[3]), 3)
        extendable_score = round(float(v[4]), 3)
        num_junctions = int(v[5])
        features.extend([bottleneck_coverage, highest_coverage, extendable_score, num_junctions])
        
        # Compatible junction coverage
        comp_junc_cov = [round(float(x) / math.log(float(sample_size)), 3) for x in v[6:6+num_junctions]]
        features.extend(calc_stats(comp_junc_cov))
        
        # Junction coverage (regardless of compatibility)
        junc_cov = [round(float(x) / math.log(float(sample_size)), 3) 
                   for x in v[6+num_junctions:6+2*num_junctions]]
        features.extend(calc_stats(junc_cov))
        
        # Cell support information
        cells_support = int(v[6+2*num_junctions])
        if cells_support == 0:
            return None
        features.append(cells_support)
        features.append(round(cells_support / float(sample_size), 3))
        
        # Cell support details
        cell_support_details = v[7+2*num_junctions:7+2*num_junctions+2*cells_support]
        cell_junction_counts = [int(cell_support_details[i+1]) 
                              for i in range(0, len(cell_support_details), 2)]
        features.extend(calc_stats(cell_junction_counts))
        
        # Number of cells covering all junctions
        cells_covering_all_junctions = sum(1 for count in cell_junction_counts 
                                         if count == num_junctions)
        features.append(round(cells_covering_all_junctions / float(sample_size), 3))
        features.append(round(cells_covering_all_junctions / float(cells_support), 3))
        
        # Fragment information
        num_fragments = int(v[7+2*num_junctions+2*cells_support])
        features.append(num_fragments)
        
        # Fragment coverage
        fragment_coverage = [round(float(x), 3) 
                           for x in v[8+2*num_junctions+2*cells_support:]]
        features.extend(calc_stats(fragment_coverage))
        
        return features
    except (IndexError, ValueError, ZeroDivisionError) as e:
        print(f"Error processing row: {e}")
        return None

def process_general_feature_file(feature_file, sample_size):
    """Process feature file and return a DataFrame."""
    feature_names = [
        'transcript_id', 'chr',
        'sample_size',
        'bottleneck_coverage', 'highest_coverage', 'extendable_score', 'num_junctions',
        'min_comp_cov', 'max_comp_cov', 'median_comp_cov', 'mean_comp_cov', 'std_comp_cov',
        'min_junc_cov', 'max_junc_cov', 'median_junc_cov', 'mean_junc_cov', 'std_junc_cov',
        'cells_support', 'ratio_cells_support',
        'min_cell_junc', 'max_cell_junc', 'median_cell_junc', 'mean_cell_junc', 'std_cell_junc',
        'ratio_cells_fl_over_sample', 'ratio_cells_fl_over_cells_support',
        'num_fragments',
        'min_frag_cov', 'max_frag_cov', 'median_frag_cov', 'mean_frag_cov', 'std_frag_cov'
    ]
    
    all_data = []
    skipped_rows = 0
    
    print(f"Processing feature file with sample size: {sample_size}")
    with open(feature_file, 'r') as f:
        reader = csv.reader(f)
        for row_num, row in enumerate(reader, 1):
            features = get_meta_feature_vector(row, sample_size)
            if features is not None:
                all_data.append(features)
            else:
                skipped_rows += 1
                if skipped_rows <= 5:
                    print(f"Warning: Skipped row {row_num} due to invalid data")
    
    if skipped_rows > 0:
        print(f"Total rows skipped: {skipped_rows}")
    
    return pd.DataFrame(all_data, columns=feature_names)

def score_general_transcripts(df, model, prob_threshold):
    """Score transcripts using the pre-trained model."""
    # Select feature columns (excluding transcript_id and chr)
    feature_cols = df.columns[2:]
    
    # Get probability scores
    probabilities = model.predict_proba(df[feature_cols])[:, 1]
    
    # Filter results based on probability threshold
    mask = probabilities >= prob_threshold
    filtered_df = df[mask].copy()
    filtered_probabilities = probabilities[mask]
    
    # Create results DataFrame with only transcripts that pass threshold
    results = pd.DataFrame({
        'transcript_id': filtered_df['transcript_id'],
        'chr': filtered_df['chr'],
        'probability_score': filtered_probabilities
    })
    
    return results

def main():
    args = parse_arguments()
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Load pre-trained model
    print(f"Loading pre-trained model from {args.model}")
    model = joblib.load(args.model)
    
    # Process input data
    print("Processing input features...")
    df = process_general_feature_file(args.input, args.sample_size)
    print(f"Processed {len(df)} transcripts")
    
    # Save processed features
    features_file = os.path.join(args.output, "general_features.csv")
    df.to_csv(features_file, index=False)
    print(f"Saved processed features to {features_file}")
    
    # Score transcripts
    print("Scoring transcripts...")
    results = score_general_transcripts(df, model, args.probability)
    
    # Save results
    output_file = os.path.join(args.output, "beaver_general_scores.csv")
    results.to_csv(output_file, index=False)
    print(f"Saved scoring results to {output_file}")
    
    # Print summary statistics
    print("\nScoring Summary:")
    print(f"Total transcripts processed: {len(df)}")
    print(f"Transcripts passing threshold: {len(results)}")
    print(f"Transcripts filtered out: {len(df) - len(results)}")
    print(f"Using probability threshold: {args.probability}")

if __name__ == "__main__":
    main()