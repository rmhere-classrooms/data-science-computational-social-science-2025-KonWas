library(igraph)
library(shiny)

# wczytanie danych
data_url <- "https://bergplace.org/share/out.radoslaw_email_email"
dfGraph <- read.csv(data_url, skip = 2, sep = " ", header = FALSE)[, 1:2]
dfGraph_original <- dfGraph
colnames(dfGraph) <- c("from", "to")

# utworzenie grafu skierowanego, simplify - pozbywa sie powtorzen
g <- graph.data.frame(dfGraph, directed = TRUE)
g <- simplify(g)

# sprawdzenie czy zgadza sie z trescia
cat("Nodes:", vcount(g), "| Edges:", ecount(g), "\n")
if (vcount(g) != 167 || ecount(g) != 5783) {
  warning("Expected 167 nodes and 5783 edges!")
}

# obliczanie wag, najpierw liczba emaili miedzy para wezlow
email_counts <- as.data.frame(table(dfGraph_original))
colnames(email_counts) <- c("from", "to", "count")
email_counts <- email_counts[email_counts$count > 0, ]

# potem suma emaili wyslanych przez kazdy wezel
total_sent <- aggregate(count ~ from, data = email_counts, FUN = sum)
colnames(total_sent) <- c("from", "total")

# obliczenie wag jako stosunek liczby emaili do sumy wyslanych
email_counts <- merge(email_counts, total_sent, by = "from")
email_counts$weight <- email_counts$count / email_counts$total

# unikalny klucz dla dopasowania wag do krawedzi grafu
email_counts$key <- paste(email_counts$from, email_counts$to, sep = "_")

# przypisanie wag do krawedzi grafu
edge_list <- as_edgelist(g, names = TRUE)
edge_keys <- paste(edge_list[,1], edge_list[,2], sep = "_")
E(g)$weight <- email_counts$weight[match(edge_keys, email_counts$key)]

# przygotowanie list sasiadow i wag krawedzi dla kazdego wezla - szybszy dostep
adj_list <- vector("list", vcount(g))
edge_weights_list <- vector("list", vcount(g))

for (i in 1:vcount(g)) {
  # wszystkie krawedzie wychodzace z wezla i
  neighbors_edges <- E(g)[.from(i)]
  if (length(neighbors_edges) > 0) {
    # sasiedzi to koncowe wezly tych krawedzi
    adj_list[[i]] <- ends(g, neighbors_edges)[, 2]
    edge_weights_list[[i]] <- E(g)$weight[neighbors_edges]
  } else {
    adj_list[[i]] <- character(0)
    edge_weights_list[[i]] <- numeric(0)
  }
}

# nazwy to np. "1" wiec zeby bylo potem latwiej sie odwolywac mapujemy je na inty
node_name_to_id <- setNames(1:vcount(g), V(g)$name)

simulate_ic <- function(g, initial_nodes, weight_mult = 1.0, max_iter = 50) {
  # zamiast korzystac z atrybutow wykorzystam listy, gdzie mamy bezposredni dostep do wartosci
  # i nie musimy przechodzi przez struktury igrapha
  # gdy wczesniej uzywalem V(g)$activated to obliczenia trwaly po kilkanascie minut
  n <- vcount(g)
  activated <- rep(FALSE, n)
  tried <- vector("list", n)
  
  activated[initial_nodes] <- TRUE
  
  # zapis historii aktywacji, ile wezlow aktywowanych w kazdej iteracji
  history <- c(length(initial_nodes))
  newly_activated <- initial_nodes
  iteration <- 1
  
  while (length(newly_activated) > 0 && iteration <= max_iter) {
    next_activated <- c()
    
    for (node_id in newly_activated) {
      neighbors <- adj_list[[node_id]]
      if (length(neighbors) == 0) next
      
      # wagi krawedzi do sasiadow oraz ich ID
      weights <- edge_weights_list[[node_id]]
      
      neighbor_ids <- node_name_to_id[neighbors]
      
      for (j in 1:length(neighbor_ids)) {
        neighbor_id <- neighbor_ids[j]
        
        # dla kazdego sasiada sprawdzenie czy nie byl juz aktywowany lub probowany przez ten wezel
        if (!activated[neighbor_id] && !(neighbor_id %in% tried[[node_id]])) {
          
          tried[[node_id]] <- c(tried[[node_id]], neighbor_id)
          
          activation_prob <- min(weights[j] * weight_mult, 1.0)
          
          # sprawdzenie czy trzeba aktywowac sasiada
          if (runif(1) < activation_prob) {
            activated[neighbor_id] <- TRUE
            next_activated <- c(next_activated, neighbor_id)
          }
        }
      }
    }
    
    # dopisanie do historii aktyacji, aktualizacja nowo aktywowanych wezlow
    history <- c(history, length(next_activated))
    newly_activated <- unique(next_activated)
    iteration <- iteration + 1
  }
  
  return(list(total = sum(activated), history = history))
}

select_nodes <- function(g, method = "random", n = NULL) {
  if (is.null(n)) n <- ceiling(vcount(g) * 0.05)
  
  if (method == "outdegree") {
    return(order(degree(g, mode = "out"), decreasing = TRUE)[1:n])
  } else if (method == "betweenness") {
    return(order(betweenness(g, directed = TRUE), decreasing = TRUE)[1:n])
  } else if (method == "closeness") {
    cls <- closeness(g, mode = "out")
    cls[is.nan(cls)] <- 0
    return(order(cls, decreasing = TRUE)[1:n])
  } else if (method == "pagerank") {
    return(order(page.rank(g, directed = TRUE)$vector, decreasing = TRUE)[1:n])
  } else {
    return(sample(1:vcount(g), n))
  }
}

cached_centralities <- list()

get_cached_nodes <- function(g, method, n = NULL) {
  if (is.null(n)) n <- ceiling(vcount(g) * 0.05)
  
  if (!(method %in% names(cached_centralities))) {
    cached_centralities[[method]] <<- select_nodes(g, method, n)
  }
  
  return(cached_centralities[[method]])
}

run_experiments <- function(g, reps = 100, weight_mult = 1.0, max_iter = 50) {
  methods <- c("outdegree", "betweenness", "closeness", "random", "pagerank")
  names <- c("Out-degree", "Betweenness", "Closeness", "Random", "PageRank")
  results <- list()
  
  # PageRank: miara centralnosci oparta na strukturze linkow - wezly sa wazne jesli
  # wskazuja na nie inne wazne wezly - osoby, ktore otrzymuja emaile od innych waznych osob
  
  # dla metod innych niz losowa mozemy obliczyc poczatkowe wezly tylko raz
  for (i in 1:(length(methods)-1)) {
    if (methods[i] != "random") {
      get_cached_nodes(g, methods[i])
    }
  }
  
  for (i in 1:length(methods)) {
    method <- methods[i]
    all_histories <- list()
    
    if (method != "random") {
      fixed_initial <- get_cached_nodes(g, method)
    } else {
      fixed_initial <- NULL
    }
    
    for (rep in 1:reps) {
      initial <- if (method == "random") select_nodes(g, method) else fixed_initial
      result <- simulate_ic(g, initial, weight_mult, max_iter)
      all_histories[[rep]] <- result$history
    }
    
    # tworzenie macierzy historii aktywacji
    max_len <- max(sapply(all_histories, length))
    history_matrix <- matrix(0, nrow = reps, ncol = max_len)
    for (rep in 1:reps) {
      hist <- all_histories[[rep]]
      history_matrix[rep, 1:length(hist)] <- hist
    }
    
    results[[method]] <- list(name = names[i], avg = colMeans(history_matrix))
  }
  
  return(results)
}


# Shinyapps
ui <- fluidPage(
  titlePanel("Information Diffusion - Independent Cascades Model"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("weight_mult", "Activation probability (% of wij):",
                  min = 10, max = 200, value = 100, step = 10, post = "%"),
      
      sliderInput("max_iter", "Max iterations:",
                  min = 1, max = 50, value = 10, step = 1),
      
      actionButton("run_sim", "Run Simulation", class = "btn-primary"),
    ),
    
    mainPanel(
      plotOutput("diffusion_plot", height = "500px"),
      br(),
      verbatimTextOutput("summary_text")
    )
  )
)

server <- function(input, output, session) {
  
  simulation_results <- eventReactive(input$run_sim, {
    run_experiments(g, 
                    reps = 100, 
                    weight_mult = input$weight_mult / 100,
                    max_iter = input$max_iter)
  })
  
  output$diffusion_plot <- renderPlot({
    results <- simulation_results()
    colors <- c("#ff0000", "#00ff00", "#0000ff", "#00f0f0", "#f000f0")
    max_len <- max(sapply(results, function(r) length(r$avg)))
    max_val <- max(sapply(results, function(r) max(r$avg)))
    
    plot(NULL, xlim = c(0, max_len - 1), ylim = c(0, max_val),
         xlab = "Iteration", ylab = "Activated nodes (average)",
         main = paste0("Information Diffusion (", input$weight_mult, "%, ", 
                       input$max_iter, " iterations)"),
         cex.main = 1.1, cex.lab = 1.0)
    
    for (i in 1:length(results)) {
      history <- results[[i]]$avg
      lines(0:(length(history) - 1), history, col = colors[i], lwd = 2.5)
    }
    
    legend("topleft", legend = sapply(results, function(r) r$name),
           col = colors, lwd = 2.5, cex = 1.0, bg = "white")
    grid()
  })
  
  output$summary_text <- renderText({
    results <- simulation_results()
    lines <- c(
      paste("Probability multiplier:", input$weight_mult, "%"),
      paste("Max iterations:", input$max_iter),
      "\nTotal activated nodes (average):"
    )
    
    for (method in names(results)) {
      total <- sum(results[[method]]$avg)
      lines <- c(lines, paste("  -", results[[method]]$name, ":", round(total, 1)))
    }
    
    paste(lines, collapse = "\n")
  })
}

shinyApp(ui = ui, server = server)