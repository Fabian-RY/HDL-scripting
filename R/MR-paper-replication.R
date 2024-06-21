#! /usr/bin/env Rscript

library(MendelianRandomization)
```{r}

```
# Generamos el input usando 2 disponibles en la propia librería
# Estos son los que se usan en el paper
MRdata_HDL_CHD <- mr_input(bx = hdlc, bxse = hdlcse, by = chdlodds, byse = chdloddsse)

# mr_ivw, mr_egger y mr_median devuelven una estimación, su error estandar e 
# intervalo de confianza dependiendo de los métodos utilizados,  junto con el 
# error residual y otros parámetros
stats <- mr_ivw(MRdata_HDL_CHD)
egger <- mr_egger(MRdata_HDL_CHD)
median <- mr_median(MRdata_HDL_CHD)

# mr_allmethods 
allmethods <- mr_allmethods(MRdata_HDL_CHD, method = "main")

mr_plot(allmethods)

mr_plot(MRdata_HDL_CHD)
