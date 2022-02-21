# bcfdelta

BCF-Delta is an application that losslessly compresses VCF/BCF files.

It differs from other "improved" VCF formats or compressors in the following ways:

  * No new format is created, instead the values inside the existing format are subtracted in a way that improves their compressability.¹
  * Works with with VCF and BCF but the improvement on BCF is much better.
  * Very simple compression scheme based on delta-compression.
  * Good at compressing files with many samples and mostly Integer-fields.
  * Low-to-Moderate memory usage (≤ 3 uncompressed records).

See the [wiki](https://github.com/h-2/bcfdelta/wiki) for a full description of the compression and a comparison with other applications.

¹ This means files can still be read with `bcftools`! Although interpretation of the shown values may depend on previous values.

## Usage

Compress a file:

```
./bcfdelta encode input_file[.vcf.gz|.bcf] output_file[.vcf.gz|.bcf]
```

Uncompress a file:

```
./bcfdelta decode input_file[.vcf.gz|.bcf] output_file[.vcf.gz|.bcf]
```

See the respective help pages (`--help`) for more details.

## Disclaimer

* This is an early preview and everything is still subject to change.
* Compression-ratio depends on input and may vary (please give us feedback!).
* No tuning has been done for (de-)compression speed, yet. It is slower than other tools (but this will improve).
* There is no support, yet, for decoding a sub-region from the file (but this will be added soon).

## Build instructions

Clone the latest main-branch:

```
cd ~/devel/                                     # or some other directory
git clone --recurse-submodules https://github.com/h-2/bcfdelta.git
```

Setup build folder:

```
mkdir -p ~/devel/bcfdelta-build/release         # or some other directory
cd ~/devel/bcfdelta-build/release
cmake -DCMAKE_BUILD_TYPE=Release ../../bcfdelta
```

Build:

```
make
```

Run:

```
./bcfdelta --help
```
