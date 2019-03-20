# Attenuation measurement files

These files are part of a project to measure material absorption coefficients.

## Absorption measurement

Includes a program to analyse transfer function or absorption coefficient from a set of recordings (`tf_from_rec.py`):

Example 

```bash
python tf_from_rec.py white_foam.aup 3 2 -f 100 4000 -e -n 32 -a 0.1
```

Use spectra calculated from white noise loops contained in `white_foam.aup` in channels **3** and **2**:

* `white_foam.aup`: audio file (only audacity for now), depends on [audacity.py](https://github.com/goiosunsw/audacity.py)
* `3 2`: use channels 3 and 2 (more channel pairs can be added for comparison)
* `-f 100 400`: excitation noise is bandlimitted to 100-4000 Hz
* `-e`: plot spectral error
* `-l 4096`: excitation noise loop is 4096 samples long
* `-n 32`: use the first 32 loops
* `-a 0.1`: microphone distance is 0.1 m
  * `-a 0.1 345`, similar but change speed of sound to 345 m/s


When comparing several measurements use:

```bash
python tf_from_rec.py audio1.aup 3 2 5 4 audio2.aup 1 0 [more options]
```

means: 
* use channels 3 and 2 from audio1.aup
* use channels 5 and 4 from audio1.aup
* use channels 1 and 0 from audio2.aup


## Transfer function measurement

```bash
python tf_from_rec.py white_foam.aup 0 1 2 -f 100 4000 -e -n 32 -a 0.1
```

Means: measure transfer function between channels 1 and 0 and between channels 2 and 0 

## Generation of excitation noise

*to be added later* 


