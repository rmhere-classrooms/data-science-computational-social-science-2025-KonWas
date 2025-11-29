library(igraph)

g <- barabasi.game(n = 1000)

summary(g)
cat("Liczba węzłów:", vcount(g), "\n")
cat("Liczba krawędzi:", ecount(g), "\n")

layout <- layout.fruchterman.reingold(g)
plot(g,
     layout = layout,
     vertex.size = 2,
     vertex.label = NA,
     edge.arrow.size = 0.2,
     main = "Graf Barabási-Albert (Fruchterman-Reingold layout)")

betweenness(g)
V(g)[betweenness(g) == max(betweenness(g))]

diameter(g)

# Różnice między grafem Barabási-Albert a Erdős-Rényi:

# Graf Barabási-Albert:
# Nowe wezly lacza sie chetniej z wezlami, ktore juz maja duzo polaczen
# Powstaja "gwiazdy" - kilka wezlow ma bardzo duzo polaczen, reszta malo
# Zblizone do prawdziwej sieci

# Graf Erdős-Rényi:
# Kazda para wezlow ma taka sama szanse na polaczenie (losowe)
# Wiekszosc wezlow ma podobna liczbe polaczen
# Model teoretyczny, w rzeczywistosci takie sieci sa rzadkie
