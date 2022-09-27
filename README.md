# Relate Python

## Motivation

[MyersGroup's Relate](https://github.com/MyersGroup/relate) is written in C++. Frankly, I am not good at C++. So, I decide to write a wrapper over this project. Maybe I will refactor them all in Python/Cython is the performance loss is acceptable.

## Build
```
$ git clone https://github.com/tcztzy/relate_py
$ cd relate_py
$ git submodule init
$ cd relate
$ cmake .
$ make
$ cd ..
$ pip install -e .
```
