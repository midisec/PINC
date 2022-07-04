# PINC
![logo](static\images\logo.png)

A powerful tool for identifying non-coding RNAs in plants by analysing k-mer frequency, cds-related features, sequence length and GC content to distinguish between the growing number of non-coding RNAs and coding RNAs in plants.



## Features

- [x] High precision (ensemble learning)
- [x] Multiple high-performance base models
- [x] Convenience of use
- [x] Automated Forecasting
- [x] Web Online

## Documents

[Documentation](http://www.pncrna.com/help)

## Get Start

> There are multiple ways to run this tool, feel free to choose one of the following method.

### Run PINC  from Web Online (Fastest)

http://www.pncrna.com/

### Run PINC  from docker (Locally、Simply)

1. Download the PINC and Add the data file to the project directory.

```bash
git clone https://github.com/midisec/PINC
cd PINC
# upload the data file (example: data.fasta)
```

> All input data must be in [fasta format](./example_data.fasta)

2. Pull and build the environment image. (Time required)

```bash
sudo docker build -t pinc_images .
```

3. Create and Enter a new container.

```bash
sudo docker run -it pinc_images bash
```

4. Execute PINC for prediction

```bash
python pinc.py -f data.fasta
```



### Run PINC from source code (Complex)

1.  Installation Environment([Autogluon](https://github.com/awslabs/autogluon)、[kentUtils](https://github.com/ENCODE-DCC/kentUtils))

2. Clone project, install related dependencies

```bash
git clone https://github.com/midisec/PINC
cd PINC
pip3 install -r requirements.txt
```


3. Execute PINC for prediction

```bash
python pinc.py -f data.fasta
```



## Usage

### Command line version

![logo](static\images\usage.png)

#### Prediction

```
python pinc.py -f <data.fasta>
```


### Website online version

#### Prediction

![logo](static\images\web_ui1.png)



After this, you will get a task page address with the uuid.

After that you can also check the history of the task by the uuid, usually it will be saved for **one month**.

![logo](static\images\web_ui2.png)

View Results and Download results

![logo](static\images\web_ui3.png)



## Contributors

<a href="https://github.com/midisec/PINC/graphs/contributors"><img src="./static/images/contributors.svg" /></a>