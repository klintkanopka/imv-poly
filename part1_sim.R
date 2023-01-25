setwd('~/projects/imv-poly')
source('imvhelper.R')
source('simhelper.R')

library(tidyverse)
library(parallel)

b_len <- 20
bs <- seq(from=-1.5, to=1.5, length.out=b_len)
N <- 100

control <- data.frame(b = c(rep(bs, each=N)))

cl <- makeCluster(20)
clusterSetRNGStream(cl, 1337)
clusterExport(cl, c('sim.fun', 'imv', 'imv_c', 'imv_t'))
sim_out <- clusterMap(cl, sim.fun,
                      b=control$b,
                      b_grid=0,
                      dgm='graded',
                      dam='graded',
                      n_i=50,
                      n_p=500,
                      SIMPLIFY=TRUE)
stopCluster(cl)


results <- prep.results(control, sim_out)

saveRDS(results, 'sim_data/sim_1_graded_graded.rds')

results |>
  select(b, starts_with('omega')) |>
  pivot_longer(-b, names_to='IMV', values_to='value') |>
  group_by(b, IMV) |>
  summarize(mean = mean(value),
            se = sd(value),
            .groups='drop') |>
  ggplot(aes(x=b, y=mean, color=IMV,
             fill=IMV)) +
  geom_vline(aes(xintercept=0), lty=2, alpha=0.5) +
  geom_line() +
  ggtitle('DGM: GRM - DAM: GRM') +
  theme_bw()

ggsave('fig/1_b_scan_graded_graded.png', width=8, height=6)




cl <- makeCluster(20)
clusterSetRNGStream(cl, 1337)
clusterExport(cl, c('sim.fun', 'imv', 'imv_c', 'imv_t'))
sim_out <- clusterMap(cl, sim.fun,
                      b=control$b,
                      b_grid=0,
                      dgm='graded',
                      dam='gpcm',
                      n_i=50,
                      n_p=500,
                      SIMPLIFY=TRUE)
stopCluster(cl)


results <- prep.results(control, sim_out)

saveRDS(results, 'sim_data/sim_1_graded_gpcm.rds')

results |>
  select(b, starts_with('omega')) |>
  pivot_longer(-b, names_to='IMV', values_to='value') |>
  group_by(b, IMV) |>
  summarize(mean = mean(value),
            se = sd(value),
            .groups='drop') |>
  ggplot(aes(x=b, y=mean, color=IMV,
             fill=IMV)) +
  geom_vline(aes(xintercept=0), lty=2, alpha=0.5) +
  geom_line() +
  ggtitle('DGM: GRM - DAM: GPCM') +
  theme_bw()

ggsave('fig/1_b_scan_graded_gpcm.png', width=8, height=6)

