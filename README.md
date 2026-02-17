## Getting Started

```bash

# Install raEDClust by source codes
git clone https://github.com/junhaiqi/raEDClust.git
cd raEDClust && make -j8  # C++11 required to compile
# Example: clustering  `test/test.fa`, and results in `test/test.clustering.out.fa`
./raEDClust test/test.fa test/test.clustering.out.fa

```


## Overview of raEDClust
raEDClust is a specialized clustering tool designed for tandem repeat units.
Unlike conventional sequence clustering tools, raEDClust employs a approximate rotation-aware edit distance metric that explicitly accounts for circular permutations of repeat units. This strategy enables more biologically meaningful clustering of tandem repeats, rotational shifts represent the same underlying structure.

Key Features:
* Designed specifically for tandem repeat unit clustering
* Uses rotation-aware edit distance to handle circular permutations
* Improves similarity estimation for tandem repeat structures
* Supports multi-threading for scalable analysis
* Direct integration with [r2rtr](https://github.com/junhaiqi/r2rtr.git) output format
## Table of contents

  * [Requirements](#requirements)
  * [Installation](#installation)
  * [Usage](#usage)
  * [Example](#example)
  * [Output](#output)
  * [Acknowledgments](#acknowledgments)
  * [License](#license)
  * [Cite](#cite)

## Requirements
raEDClust runs Linux and requires gcc 9.3.0+.


## Installation

```bash
git clone https://github.com/junhaiqi/raEDClust.git
cd raEDClust
make -j8  # C++11 required to compile
```
Then, there will be a binary file called raEDClust.

## Usage

Basic command:
```bash
rotCluster <in.fa> <out.fa> [Options:]
```
The specific parameters are as follows:
```bash
Options: 
  -l FLOAT   Minimum length ratio (default: 0.85)
  -s FLOAT   Similarity threshold (default: 0.90)
  -c STR     Output clustering results to this file (default: None)
  -x BOOL    in.fa is in r2rtr format (default: 0)
  -t INT     Number of threads (default: max available)
  -h         Show this help message
```

## Example

Clustering  `test/test.fa`, and results in `test/test.clustering.out.fa`:

```bash
./raEDClust test/test.fa test/test.clustering.out.fa
```

## Output
The output is a FASTA file containing the clustered tandem repeat units. Each record follows the format below:

```bash
>cluster1|S1_2!>chrX:57819763-60927195!1334283!1360573!-:2057:2119:64
TTTTCGCAGAATCTGCAAGTGGACATTTGGAGCGCTTTCAGGCCTGTGGTGGAAAAGGCCTGAAAGCCTT...
```
The FASTA header contains the following fields:
```bash
>ClusterID!OriginalReadID!GenomicRegion!Start!End!Strand:InferredLength:ConsensusLength:SupportCount
```
Using the example above:

* cluster1
Cluster ID.

* S1_1!>chrX:57819763-60927195!2042180!2063944!+
r2rtr format: The original read identifier together with its mapped genomic region, including chromosome coordinates and strand orientation.

* 2057
r2rtr format: Inferred tandem repeat unit length derived from read-to-read overlap distance patterns.

* 2119
r2rtr format: Final consensus unit length after sequence refinement.

* 64    
r2rtr format: Number of repeat fragments supporting the consensus unit.


## License 
MIT License.

## Cite
None

