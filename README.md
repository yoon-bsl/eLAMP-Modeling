# eLAMP-Modeling
Mathematical modeling of the creation and adsorption of eLAMP amplicons to the oil-water interface  

Outputs a figure depicting the percentage of emulsified droplet surface area saturated with LAMP amplicons over time.

<i>Based on the model found in Ulep, T.-H., Day, A.S., Sosnowski, K., Shumaker, A., & Yoon, J.-Y. Interfacial Effect-Based Quantification of Droplet Isothermal Nucleic Acid Amplification for Bacterial Infection. Sci. Rep. 9, 9629 (2019).</i>

---
Example command line run:  
```python
python main.py -l 193 -t med -e exp -d 22 -a 10
```

Input arguments:
* l - the base pair length of the target gene
* t - type of statistic used for droplet measurement (med for median droplet size, or avg for average)
* e - exponential (exp) or logarithmic (log) amplicon creation model
* d - doubling time (in seconds) of the LAMP reaction
* a - number of different amplicon sizes to investigate (growing geometrically from the target base pair length) 

In order to run on a given machine, clone this repository to a directory of your choice, and use pip installer to download the necessary Python libraries used in this work:  
```
pip install -r requirements.txt
```

Example model output given the above command line instance:  

<img align="left" src="src\ExampleOutput.png" width="500px">
