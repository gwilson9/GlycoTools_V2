<h1 align="center">
  <br>
    <a><img src="https://i.ibb.co/G2ZfVkH/20191213-Glyco-Tools-11.png" alt="20191213-Glyco-Tools-11" border="0" width = "1000"></a>
  <br>
</h1>


# GlycoTools

GlycoTools simplifies the processing of glycoproteomic data generated by the Byonic™ search engine.

### Key Features

* Organization of GlycoPSMs into Glycopeptide, Glycosite, Glycoprotein, and glycan data tables for user defined conditions and replicates. 
* Filter glycoPSMs by user-defined thresholds.
* Protein inference extended across all data files in an experiment.
* Glycosite localization for all N- and O-sites on a peptide backbone.
* Detection of GlycoPSMs that are likely artifacts of in-source fragmentation of a larger glycopeptide precursor.
* Visualization of peptide elution, annotated spectra, and other relevant metrics of data quality.
* Intuitive user experience in a contained GUI.

### Prerequisites

GlycoTools is developed in C# with Microsoft Visual Studio 2015 and the Microsoft.NET Framework version 4.5.2 and should be amenable with Windows 7 or later.

### Usage

Drag and drop functionality allows the user to add Byonic result (.byrslt) and Thermo RAW files into the 'Data Upload' table. Raw files will be paired with the corresonding result file by matching file names, which are conserved by Byonic convention. The user then defines the condition and relicate affiliations of each file in the experiment. The user can define quality filters and optional features, including protein inference, glycan localization, and detection of in-source fragments, before running the program. A progress log will provide updates as GlycoTools processes the data. Once finished, the tables and charts in subsequent tabs will be populated with the processed data. Processed results can be dumped into .csv files from the 'Result Tables' tab or exported manually from the associated SQLite database located in the user-defined output folder. This database can also be dropped into a fresh instance of GlycoTools to view previously processed data.

### Installation


