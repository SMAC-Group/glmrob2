## Binary response --------------
data(vaso)

## classical fit
Vfit1 <- glm(Y ~ log(Volume) + log(Rate), family=binomial, data=vaso)
coef(Vfit1)

## robust fit (robustbase)
Vfit2 <- glmrob(Y ~ log(Volume) + log(Rate), family=binomial, data=vaso,
                method="Mqle", control=glmrobMqle.control(tcc=3.5))
coef(Vfit2)

## robust fit 2 (new method, allows for different weight functions)
Vfit3 <- glmrob2(Y ~ log(Volume) + log(Rate), family=binomial, data=vaso)
coef(Vfit3)
