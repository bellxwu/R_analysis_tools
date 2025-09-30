# Description: Make boxplots for comparisons between two groups w/ or w/o facets
# Authour: Bell Wu
# Date created: 2025.09.25
# df = long df created from pivot_longer

make_boxplot = function(df, x, y, group = NULL, facet) {
  box = ggplot(df, mapping = aes(x = as.name(x), y = as.name(y), fill = as.name(group))) +
    geom_boxplot() +
    theme_few() +
    stat_compare_means(label = "p.format",
                       vjust = 1,
                       bracket.size = ) +
    facet_wrap(facet)
  return(box)
}

df = df_long
x = "Sample"
y = "TPM"
fill = "Sample"
facet = "Gene"

make_boxplot(df_long, x, y, group = fill, facet)
