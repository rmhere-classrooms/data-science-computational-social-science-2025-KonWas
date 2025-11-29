library(igraph)

g <- erdos.renyi.game(n = 100, p.or.m = 0.05)

summary(g)
# Graf nie jest ważony bo w podsumowaniu nie ma litery 'W'

V(g)
E(g)

E(g)$weight <- runif(length(E(g)), 0.01, 1)

summary(g)
# Graf jest ważony bo w podsumowaniu jest litera 'W'

node_degrees <- degree(g)
node_degrees

hist(node_degrees,
     main = "Rozkład stopni węzłów w grafie Erdős-Rényi",
     xlab = "Stopień węzła",
     ylab = "Liczba węzłów",
     col = "lightblue",
     breaks = 10)

cl <- components(g)
cat("Liczba klastrów:", cl$no, "\n")
cat("Rozmiary klastrów:", cl$csize, "\n")

pr <- page_rank(g)$vector

plot(g,
     vertex.size = pr * 500,
     vertex.label = NA,
     edge.arrow.size = 0.3,
     main = "Graf Erdős-Rényi z rozmiarami węzłów wg PageRank")
