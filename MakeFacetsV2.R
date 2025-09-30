# Title: Facet bar plots for marker analyses
# Description: Generates a facet plots comparing two different variables
# Date created: 2025.09.26

# input_df = input data frame
# input_cols = variables that will be analyzed into an individual facet plots (ie. markers to be analyzed)
# long = if TRUE, means the input df is already in a long df format
# x = independent variable that will be plotted on bar plots, can be either "Treatment" or "Replicates"
# y = dependent variable that will be plotted on bar plots
# facet_var = How the plots will be faceted from 
# plot_name = Title of plot

# Changelog: (from version 1)
# 2025.09.26 - Added in argument where you can choose the plot
# Added in argument where you can choose a variable to facet
# Added a base ggthemes theme
# Added a way to check if the df is a long df or short df
# 2025.09.27 - Add optional stat-test entry
# Added 

MakeFacets <- function(input_df,
                       long = FALSE,
                       input_cols, 
                       x = "x", # labels for x
                       y = "y", # labels for y
                       fill, # variable you want to add colour to
                       plot_type,
                       intervals = 11,
                       add_stat = TRUE,
                       stat_test = "t.test",
                       facet_var = "None") {
# ------------------------------------------------------------------------------  
# Option whether input df is long format, if not will convert into long
  if (!long) {
    df_long = input_df %>%
      pivot_longer(cols = input_cols,
                   names_to = x,
                   values_to = y)
  } else {
    df_long = input_df
  }
# ------------------------------------------------------------------------------
# Creating the base data needed for the ggplot: 
  plots = ggplot(df_long, mapping = aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill]]))
# ------------------------------------------------------------------------------ 
# Bar plot format
  bar = stat_summary(fun = median,
                     geom = "bar",
                     width = 0.6)
# Box plot format
  box = geom_boxplot()
# Histogram format
  hist = geom_bar(position = "identity")
# ------------------------------------------------------------------------------
# Choice of either bar, box, or histogram: 
  if (plot_type == "bar") {
  plots = plots + bar
  } else if (plot_type == "box") {
    plots = plots + box
  } else if (plot_type == "hist") {
    if (any(colnames(df_long) == "bin")) {
      plots = ggplot(df_long, mapping = aes(x = bin, fill = .data[[fill]])) +
        hist
    } else {
      breaks = seq(0, max(df_long$TPM) + 1, length.out = intervals)  # create bins
      df_long = df_long |> 
        mutate(bin = cut(df_long$TPM, breaks = breaks, right = FALSE)) # add bins to df_long
      plots = ggplot(df_long, mapping = aes(x = bin, fill = .data[[fill]])) +
        hist
    }
  } else {
    stop("Please enter plot type: bar, box, hist")
  }
# ------------------------------------------------------------------------------ 
# add optional statistical test
if (add_stat) {
  plots = plots +
    stat_compare_means(data = df_long,
                       method = stat_test,
                       label = "p.format",
                       size = 4,
                       vjust = 1)
} else {
  plots = plots
}
# ------------------------------------------------------------------------------   
# add the faceting as an optional variable
  if (facet_var == "None") {
    plots = plots
  } else {
    plots = plots + facet_wrap(~ .data[[facet_var]],
                               scales = "free_y")
  }
# ------------------------------------------------------------------------------  
# added in theme  
  plots = plots + 
    theme_few() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))
  return(plots)
}


### Random testing ------

