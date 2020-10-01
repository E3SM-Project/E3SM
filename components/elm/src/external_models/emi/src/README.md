### Instructions to run the demo code

#### 

1. Specify the Fortran compiler 
```
export FC=<fortran-compiler>
```

2. Run configure command that will create a new directory  `./build/`
```
make FC=$FC config
```

3. Install the demo code that will create a new directory `./local`.
```
make FC=$FC install
```

4. Run the demo code
```
cd local/bin
./demo
```
