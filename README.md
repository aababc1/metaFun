# metaFun : pipeline for metagenomic data analysis
[![DOI](https://img.shields.io/badge/DOI-10.1080%2F19490976.2025.2611544-blue)](https://doi.org/10.1080/19490976.2025.2611544)
[![Conda](https://img.shields.io/conda/v/bioconda/metafun)](https://anaconda.org/bioconda/metafun)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Documentation](https://img.shields.io/badge/docs-readthedocs-blue)](https://metafun-doc.readthedocs.io/)

A scalable and agile analysis pipeline for **meta**genomic big data with fast and unified **Fun**ctional searches.


<img width="217" height="195" alt="image" src="https://github.com/user-attachments/assets/ad33f6ef-6ee9-477b-8fc2-bd18ce205baf" />

<img width="2550" height="3300" alt="pipeline_flowchart renew" src="https://github.com/user-attachments/assets/1e673042-d22b-4574-a511-c98ac5df8e38" />




# Quick Start

## Installation
```bash
conda create -c bioconda -c conda-forge -n metafun bioconda::metafun
conda activate metafun
```

### If you have mamba, you can install metaFun with mamba
mamba create -c bioconda -c conda-forge -n metafun metafun 

### Download databases
```bash
metafun -module DOWNLOAD_DB
```

### Run analysis
```bash
metafun -module RAWREAD_QC  -i input_reads/
```




# Documentation 📖 
https://metafun-doc.readthedocs.io/en/latest/index.html



# FAQ  
Please  use "Issues tab in this page"  or  ask in slack .
https://join.slack.com/t/metafun1/shared_invite/zt-3jud1l51y-fyhyco0tuRlxODet~O6~7w


## Support & Contact

metaFun is actively maintained with ongoing feature development and bug fixes.
There is no restriction to use this program.

If you have any questions, feature requests, or encounter issues, feel free to reach out:
- **GitHub Issues**: https://github.com/aababc1/metaFun/issues
- **Email**: aababc1@yonsei.ac.kr  
- **Slack**: [Join our workspace](https://join.slack.com/t/metafun1/shared_invite/zt-3jud1l51y-fyhyco0tuRlxODet~O6~7w)

When reporting bugs, please include your OS, the module/step where the error occurred, and relevant log files (`.nextflow.log`, `.command.log`).


## 📚 Citation

If you use metaFun in your research, please cite:

> Lee HG, et al. (2026). metaFun: An analysis pipeline for metagenomic big data with fast and unified functional searches. *Gut Microbes*, 17(1). https://doi.org/10.1080/19490976.2025.2611544

[![DOI](https://img.shields.io/badge/DOI-10.1080%2F19490976.2025.2611544-blue)](https://doi.org/10.1080/19490976.2025.2611544)

<details>
<summary>BibTeX</summary>
```bibtex
@article{lee2026metafun,
  title={metaFun: An analysis pipeline for metagenomic big data with fast and unified functional searches},
  author={Lee, Hyeon Gwon and others},
  journal={Gut Microbes},
  volume={17},
  number={1},
  year={2026},
  publisher={Taylor \& Francis},
  doi={10.1080/19490976.2025.2611544}
}
```
</details>
