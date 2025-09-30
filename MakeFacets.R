# Title: Facet bar plots for marker analyses
# Description: Generates a facet plots comparing two different variables
# Date updated: 08.06.2025

# input_df = input data frame
# input_cols = variables that will be analyzed into an individual facet plots (ie. markers to be analyzed)
# x_name = name of the x variable plotted; the independent converted from the input columns
# y_name = name of the y variable plotted; the dependent value of the input x_name
# x = independent variable that will be plotted on bar plots, can be either "Treatment" or "Replicates"
# y = dependent variable that will be plotted on bar plots
# facet_var = How the plots will be faceted from 
# plot_name = Title of plot


MakeFacets <- function(input_df, input_cols, 
                       x_name = "Marker", y_name = "Proportion", 
                       x = "Marker", 
                       y = "Proportion",
                       facet_var, stat_test = "t.test",
                       plot_name = "Title") {
# Create the long format of the input_df  
  df_long <- input_df %>%
    pivot_longer(
      cols = input_cols,
      names_to = x_name,
      values_to = y_name
    )
# Create single bar plot  
  plots <- ggplot(df_long, mapping = aes_string(x = as.name(x), y = as.name(y))) +
    stat_summary(
      fun = median,
      geom = "bar",
      width = 0.6,
      fill = "grey",
      color = "black"
      ) +
  geom_point(
    size = 1,
    color = "black",
    alpha = 0.8
      )
  
# Create facet plots
  facets <- plots + 
    facet_wrap(~.data[[facet_var]],
               nrow = 3,
               scales = "free_y"
    ) + 
    # adjust the y axis to make it bigger for each facet plot
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.35))
    ) +
    # add the titles and theme for the aesthetics  
    labs(title = plot_name) +
    theme_linedraw(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
# Add statistical test
  calc_plots <- facets +
    stat_compare_means(
      method = stat_test,
      label = "p.format",
      size = 4,
      vjust = -1.5
    )
  return(calc_plots)
}


### Random testing ------
# HuH82_S_long <- H82_S_prop %>%
#   pivot_longer(
#     cols = all_of(cols),
#     names_to = "Marker",
#     values_to = "Proportion"
#   )
# 
# S_plots <- ggplot(HuH82_S_long, mapping = aes(x = Treatment, y = Proportion)) +
#   stat_summary(
#     fun = median,
#     geom = "bar",
#     width = 0.6,
#     fill = "grey",
#     color = "black"
#   ) +
#   geom_point(
#     size = 1,
#     color = "black",
#     alpha = 0.8
#   ) +
#   scale_y_continuous(
#     expand = expansion(mult =  c(0, 0.1))
#   )
# print(S_plots)
# 
# S_facets <- S_plots +
#   facet_wrap(~Marker,
#              nrow = 3,
#              scales = "free_y"
#   ) + 
#   # adjust the y axis to make it bigger for each facet plot
#   scale_y_continuous(
#     expand = expansion(mult = c(0, 0.35))
#   ) +
#   # add the titles and theme for the aesthetics  
#   labs(title = "ATMi+IR vs ATMi+IR+aPDL1") +
#   theme_linedraw(base_size = 13) +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
#   )
# print(S_facets)
# 
# calc_plots <- S_facets +
#   stat_compare_means(
#     method = "t.test",
#     label = "p.format",
#     size = 4,
#     vjust = -1.5
#   )
# print(calc_plots)

# 
# input_df <- H82_2comb
# input_cols <- cols
# x_name <- "Marker"
# y_name <- "Proportion"
# x <- "Sample"
# y <- "Proportion"
# facet_var <- "Marker"
# plot_name = "ATMi+IR: Spleen vs TINF"