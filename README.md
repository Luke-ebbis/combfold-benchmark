# Combfold benchmark

The pipeline downloads selected structures from the PDB, and produces their
Combfold predictions. These Combfold predictions are then aligned to
the known structure and the USalign TM score is reported. The Combfold
structures can be found in `results/data/`, and the TM scores in
`results/align`.

## Installation

The following commands install dependencies of the workflow in the _current_
directory within the `.pixi` folder. After installation, you cannot move the
folder without re-installling all the dependencies. If there is a problem with a
corrupted environment, you should remove the `.pixi` and the `.snakemake`
folders.

```bash
curl -fsSL https://pixi.sh/install.sh | bash
# ... cd <this repo>
pixi install
```

## Usage

Run `pixi run help` for the help page.

The proteins that are part of the benchmark are listed in `config/config.yml`:

```toml
colabfold:
  number_of_models: "5"
assess:
  - "6dv2-assembly1"
  - "6K71"
  - "5TD9"
  - "5WBJ"
```

Under the `assess`, it is possible to specify IDs from RCSB PDB. If there
are multiple assemblies, one can specify assembly#.

## Details

![](resources/pipeline.png)
