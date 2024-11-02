# Introduction

Beaver implements an efficient algorithm to provide more accurate transcript assembly at single-cell resolution.
The datasets and scripts used to compare the performance of Beaver with other assemblers are available at
[beaver-test](https://github.com/Shao-Group/beaver-test).

# Installation

Clone the latest Beaver repository:

```bash
git clone https://github.com/Shao-Group/beaver.git
cd beaver
```

Generate configuration files:

```bash
aclocal        # Generate aclocal.m4
autoconf       # Generate configure script
autoheader     # Generate config.h.in
automake -a    # Generate Makefile.in
```

Configure and compile:

```bash
./configure    # Generate Makefile
make          # Compile the source code
```

The executable will be generated at `/path/to/your/beaver/src/beaver`.

# Usage

```
./beaver -i <input-gtf-list> -o <output-prefix> [options]
```

A quick start example is available at `beaver/sample_data/`. To play with this example:

```bash
./beaver ../sample_data/sample_input.list ../sample_data/beaver_test > test.log
```

Output will be generated at `./sample_data/` with prefix `beaver_test`.

## Format of Input and Output

### Input Individual Assembly List

 Each line of `input-gtf-list` points a path to the individual assembly (GTF format) of a single cell or sample; see [sample_input.list](https://github.com/Shao-Group/beaver/blob/master/sample_data/sample_input.list).

### Output Files Structure

1. **Cell-Specific Assemblies**
   * Location: `output-prefix_sgtf/X.gtf`
   * X ranges from 1 to number of cells
2. **Meta Assembly**
   - Location: `output-prefix.gtf`
   - Contains collective transcripts from all cells
3. **Feature Files**
   - General features (`output-prefix_feature.csv`)
   - Cell-specific features (`output-prefix_sgtf/X_feature.csv`); X ranges from 1 to number of cells

## Options

| Parameters | Type    | Default Value | Description                              |
| ---------- | ------- | ------------- | ---------------------------------------- |
| -t         | integer | 1             | Number of threads.                       |
| -c         | double  | 0             | Minimum transcript output coverage (0-1) |
| -pn        | integer | 15            | Maximum number of paths per node.        |
| -pg        | integer | 100           | Maximum number of paths per graph.       |

# Scoring Cell-Specific Transcripts with Random Forest Model

Beaver employs two-stage random forest pipeline (`Beaver_General` and `Beaver_Specific`) for scoring transcripts, available for download from [Zenodo](https://doi.org/10.5281/zenodo.14014750). Use the provided Python script `beaver_general.py` and `beaver_specific.py` with the models. Please note that if only care about the global expression of cell populations in the dataset, `Beaver_General` can work alone without `Beaver_Specific`.

## Dependencies

Required Python libraries: numPy, pandas, scikit-learn, joblib

- Using pip:

  ```bash
  pip install numpy pandas scikit-learn joblib
  ```

- Using conda (recommended):

  ```bash
  conda install numpy pandas scikit-learn joblib
  ```

## Usage

Score transcripts with the syntax below:

```bash
python3 Beaver_General.py -s -i <general_transcript_feature_csv> -m <pretrained_general_model.joblib> -p <min_probability_score> -o <output_general_dir>
```

```bash
python3 Beaver_Specific.py -i <specific_transcript_feature_dir> -g <output_general_dir> -m <pretrained_specific_model.joblib> -p <min_probability_score> -o <output_specific_dir>
```

| Parameter | Type   | Default | Description                                                  |
| --------- | ------ | ------- | ------------------------------------------------------------ |
| -i        | String |         | Input general feature CSV file or input specific feature directory |
| -m        | String |         | Path to the pre-trained model file for scoring.              |
| -n        | Integer |         | Number of cells/samples.              |
| -o        | String |         | Output directory                                             |
| -p        | String | 0.2     | Minimum probability score threshold (range: 0 to 1).         |
| -s        |        |         | Save processed feature files and estimated scores in the output directory. |
| -g        | String |         | Directory of  Beaver-General's results; feed into Beaver-Specific. |

## Example

Here is an example to run the scoring pipeline:

```
python3 beaver_general.py -s -n 2 -i sample_data/beaver_test_feature.csv -m model_real_HEK293T_meta_roc=0.824.joblib -p 0.2 -o sample_data/beaver_test_general
```

```
python3 beaver_specific.py -s -n 2 -i sample_data/beaver_test_sgtf -m model_real_HEK293T_specific_roc=0.796.joblib -p 0.6 -o sample_data/beaver_test_specific
```

