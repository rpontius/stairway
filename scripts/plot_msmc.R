suppressMessages({
  library(tidyverse)
})


pid = snakemake@params[[1]]
generations = snakemake@params[[3]]
mutation_rate = snakemake@params[[2]]
msmc_input = snakemake@input[[1]]
msmc_input_split = snakemake@input[[2]]
plot_output = snakemake@output[[1]]
plot_output_split = snakemake@output[[2]]

data<-read.table(msmc_input, header=TRUE)
plot.title = pid
data %>%
  #filter(time_index <= 57) %>%
  mutate(max_time = max(time_index - 6)) %>% 
  filter(time_index <= max_time) %>%
  mutate(YearsAgo = (left_time_boundary/mutation_rate)*generations,
         Ne = (1/lambda)/(2*mutation_rate)) %>%
  ggplot(aes(x = YearsAgo, y = Ne)) +
  	geom_step() +
  	scale_x_log10() +
  	theme_bw() + ggtitle(plot.title)
ggsave(plot_output, width = 8, height = 6)

data<-read.table(msmc_input_split, header=TRUE)
plot.title = pid
data %>%
  #filter(time_index <= 57) %>%
  mutate(max_time = max(time_index - 6)) %>% 
  filter(time_index <= max_time) %>%
  mutate(YearsAgo = (left_time_boundary/mutation_rate)*generations,
         Ne = (1/lambda)/(2*mutation_rate)) %>%
  ggplot(aes(x = YearsAgo, y = Ne)) +
  	geom_step() +
  	scale_x_log10() +
  	theme_bw() + ggtitle(plot.title)
ggsave(plot_output_split, width = 8, height = 6)