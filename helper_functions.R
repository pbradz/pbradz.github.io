### Helper functions


min_nonzero <- function(x) min(x[x > 0])

ord_plot_wrapper <- function(ord,
                             attribute,
                             title="ordination",
                             cols=palette()[c(3, 2, 4)]) {
  ord_colors <- cols[attribute]
  plot(ord, display = "sites", type = "n", main = title)
  points(ord, display = "sites", col = ord_colors, bg = ord_colors, pch = 21)
  ordiellipse(ord, attribute, col = cols, lwd = 3, draw = "polygon", label = TRUE)
}

reorder_rows <- function(tbl, col, match_to) {
  tbl[match(tbl[[col]], match_to), ]
}

# For some reason this isn't a function in tibble
as_data_matrix <- function(tbl, j=1) {
  m <- data.matrix(tbl %>% select(-!!(enquo(j))))
  rownames(m) <- deframe(tbl %>% select(!!(enquo(j))))
  m
}

nw <- function(x) names(which(x))

# Let's define a function to do the Fisher test. We provide a "background" set
# which is the set of genes tested.
do_fisher <- function(list1, list2, background) {
  l1 <- intersect(list1, background)
  l2 <- intersect(list2, background)
  both <- intersect(l1, l2) %>% length
  only1 <- setdiff(l1, l2) %>% length
  only2 <- setdiff(l2, l1) %>% length
  neither <- setdiff(background, union(l1, l2)) %>% length
  cont.table <- matrix(nr=2, c(both, only1, only2, neither))
  return(fisher.test(cont.table, alt = "g"))
}


# This is a wrapper to do the actual enrichment, given a mapping of genes to
# pathways and a set of hits.
enrich <- function(hits, mapping, background, minsize = 1) {
    m.hits <- intersect(hits, background)
    mapping.reduced <- filter(mapping, .data[[colnames(mapping)[1]]] %in% background)
    possible <- mapping.reduced[, 1] %>% deframe
    pways.total <- mapping.reduced[, 2] %>% distinct %>% deframe
    pw.lengths <- sapply(pways.total, function(y) {
        (mapping[,2] == y) %>% which %>% (function(z) {
            intersect(mapping[z, 1] %>% deframe, background) }) %>%
            length })
    pways <- nw(pw.lengths >= minsize)
    enrichments <- vector("numeric", length(pways))
    overlaps <- list()
    hits.total <- length(m.hits)
    sets <- list()
    for (np in 1:length(pways)) {
        p <- pways[np]
        names(enrichments)[np] <- p
        in.pway <- intersect(background, mapping[which(mapping[, 2] == p), 1] %>% deframe)
        ft <- do_fisher(in.pway, m.hits, background)
        pway.nohits <- length(setdiff(in.pway, m.hits))
        pway.hits <- length(intersect(in.pway, m.hits))
        nopway.hits <- length(setdiff(m.hits, in.pway))
        nopway.nohits <- length(setdiff(background, union(m.hits, in.pway)))
        cont.table <- (matrix(nr=2, c(pway.hits, nopway.hits, pway.nohits,
                                      nopway.nohits)))
        enrichments[np] <- ft$p.value
        overlaps[[p]] <- intersect(in.pway, m.hits)
        sets[[p]] <- list(PH = pway.hits, PNH = pway.nohits, NPH = nopway.hits,
                          NPNH = nopway.nohits)
    }
    return(list(enrichments = enrichments, overlaps = overlaps, sets = sets))
}

annotate_enr <- function(enr, descs, cutoff = 0.25) {
    qvals <- qvalue(enr)$qvalues
    sigs <- qvals[qvals <= cutoff]
    sig_tbl <- enframe(sigs, name=colnames(descs)[1], value="qval") %>%
        arrange(desc(qval)) %>%
        mutate(qval=format(qval, digits=3, scientific=TRUE))
    results <- left_join(sig_tbl,
                         descs)
    results
}
