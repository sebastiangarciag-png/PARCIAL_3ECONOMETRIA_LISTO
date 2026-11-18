install.packages("WDI")
library(WDI)
datos <- WDI(country = "CO", start = 1990, end = 2024,
             indicator = c(consumo = "NE.CON.PRVT.KD", pib = "NY.GDP.MKTP.KD",
                           tasa_real = "FR.INR.RINR", inflacion = "FP.CPI.TOTL.ZG",
                           desempleo = "SL.UEM.TOTL.ZS"))
write.csv(datos, "datos_colombia.csv", row.names = FALSE)
head(datos)
summary(datos)
summary(datos)

options(timeout = 300)
library(WDI)

cons <- WDI(country = "CO", start = 1990, end = 2024,
            indicator = c(consumo = "NE.CON.PRVT.KD"))
tasa <- WDI(country = "CO", start = 1990, end = 2024,
            indicator = c(tasa_real = "FR.INR.RINR"))

base <- merge(datos[, c("year","pib","inflacion","desempleo")],
              cons[, c("year","consumo")], by = "year")
base <- merge(base, tasa[, c("year","tasa_real")], by = "year")
summary(base)

write.csv(base, "base_consumo_colombia.csv", row.names = FALSE)

# ============================================================
# Parcial III - Econometría I
# Función de consumo de Colombia, 1991-2024 (datos: Banco Mundial, WDI)
# Modelo: ln(consumo) = b0 + b1*ln(PIB) + b2*tasa_real + b3*inflacion + b4*desempleo + u
# ============================================================
# Paquetes (si falta alguno: install.packages(c("lmtest","car","tseries")))
library(lmtest)    # bptest, dwtest, bgtest, resettest
library(car)       # vif
library(tseries)   # jarque.bera.test

dir.create("figuras", showWarnings = FALSE)
dir.create("tablas",  showWarnings = FALSE)
sink("resultados.txt", split = TRUE)   # guarda todo lo que se imprime en resultados.txt

# ---------- 0. Datos ----------
base <- read.csv("base_consumo_colombia.csv")
base <- na.omit(base[, c("year", "consumo", "pib", "tasa_real", "inflacion", "desempleo")])
base$lconsumo <- log(base$consumo)
base$lpib     <- log(base$pib)
cat("\nObservaciones usadas:", nrow(base), " (", min(base$year), "-", max(base$year), ")\n")

vars <- c("lconsumo", "lpib", "tasa_real", "inflacion", "desempleo")
X <- base[, vars]

# ---------- (b) Análisis descriptivo ----------
cat("\n===== (b) ESTADÍSTICOS DESCRIPTIVOS =====\n")
desc <- data.frame(
  media   = sapply(X, mean),
  mediana = sapply(X, median),
  desv    = sapply(X, sd),
  minimo  = sapply(X, min),
  maximo  = sapply(X, max),
  cv      = sapply(X, function(z) sd(z) / abs(mean(z)))   # coef. de variación
)
print(round(desc, 3))
write.csv(round(desc, 3), "tablas/descriptivos.csv")

cat("\n----- Matriz de correlaciones -----\n")
cors <- round(cor(X), 3)
print(cors)
write.csv(cors, "tablas/correlaciones.csv")

# Gráficos
png("figuras/01_series_tiempo.png", width = 1400, height = 1000, res = 150)
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
for (v in vars) plot(base$year, base[[v]], type = "l", lwd = 2, col = "steelblue",
                     xlab = "Año", ylab = v, main = paste("Evolución de", v))
dev.off()

png("figuras/02_histogramas.png", width = 1400, height = 1000, res = 150)
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
for (v in vars) hist(base[[v]], col = "lightblue", border = "white",
                     main = paste("Histograma de", v), xlab = v)
dev.off()

png("figuras/03_cajas.png", width = 1400, height = 1000, res = 150)
par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
for (v in vars) boxplot(base[[v]], col = "lightblue", main = paste("Caja de", v), ylab = v)
dev.off()

png("figuras/04_dispersion_vs_consumo.png", width = 1400, height = 1000, res = 150)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
for (v in c("lpib", "tasa_real", "inflacion", "desempleo")) {
  plot(base[[v]], base$lconsumo, pch = 19, col = "steelblue",
       xlab = v, ylab = "lconsumo", main = paste("lconsumo vs", v))
  abline(lm(base$lconsumo ~ base[[v]]), col = "red", lwd = 2)
}
dev.off()

png("figuras/05_matriz_dispersion.png", width = 1400, height = 1400, res = 150)
pairs(X, pch = 19, col = "steelblue", main = "Matriz de diagramas de dispersión")
dev.off()

# ---------- (c) Estimación por MCO ----------
cat("\n===== (c) ESTIMACIÓN MCO =====\n")
modelo <- lm(lconsumo ~ lpib + tasa_real + inflacion + desempleo, data = base)
print(summary(modelo))
cat("\nIntervalos de confianza al 95%:\n")
print(round(confint(modelo), 4))
write.csv(round(summary(modelo)$coefficients, 5), "tablas/coeficientes_MCO.csv")

# ---------- (d) Contrastes de especificación ----------
cat("\n===== (d) CONTRASTES DE SUPUESTOS DE GAUSS-MARKOV =====\n")

cat("\n--- Linealidad / especificación: RESET de Ramsey (potencias 2 y 3 del ajustado) ---\n")
print(resettest(modelo, power = 2:3, type = "fitted"))

cat("\n--- Multicolinealidad: VIF (regla: VIF > 10 problema grave; > 5 atención) ---\n")
print(round(vif(modelo), 3))

cat("\n--- Heterocedasticidad: Breusch-Pagan ---\n")
print(bptest(modelo))

cat("\n--- Heterocedasticidad: White (versión con ajustado y ajustado^2) ---\n")
print(bptest(modelo, ~ fitted(modelo) + I(fitted(modelo)^2), data = base))

cat("\n--- Autocorrelación: Durbin-Watson ---\n")
print(dwtest(modelo))

cat("\n--- Autocorrelación: Breusch-Godfrey (orden 1 y orden 2) ---\n")
print(bgtest(modelo, order = 1))
print(bgtest(modelo, order = 2))

cat("\n--- Normalidad de residuos: Jarque-Bera ---\n")
print(jarque.bera.test(residuals(modelo)))

# Gráficos de diagnóstico
png("figuras/06_diagnostico_residuos.png", width = 1400, height = 1000, res = 150)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
plot(fitted(modelo), residuals(modelo), pch = 19, col = "steelblue",
     xlab = "Valores ajustados", ylab = "Residuos", main = "Residuos vs ajustados")
abline(h = 0, col = "red", lty = 2)
qqnorm(residuals(modelo), pch = 19, col = "steelblue", main = "Q-Q de residuos"); qqline(residuals(modelo), col = "red")
hist(residuals(modelo), col = "lightblue", border = "white", main = "Histograma de residuos", xlab = "Residuos")
acf(residuals(modelo), main = "ACF de residuos")
dev.off()

sink()
cat("\nListo. Resultados en resultados.txt, gráficos en /figuras y tablas en /tablas\n")