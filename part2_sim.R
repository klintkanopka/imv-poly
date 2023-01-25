setwd('~/projects/imv-poly')
source('imvhelper.R')
source('simhelper.R')

library(tidyverse)
library(parallel)

n_cores <- 20

b_len <- 40
bs <- seq(from=-2, to=2, length.out=b_len)
N <- 100

control <- data.frame(b = c(rep(bs, each=N)))

cl <- makeCluster(n_cores)
clusterSetRNGStream(cl, 1337)
clusterExport(cl, c('sim.fun', 'imv', 'imv_c', 'imv_t'))
sim_out <- clusterMap(cl, sim.fun,
                      b=control$b,
                      b_grid=c(-1,0,1),
                      dgm='graded',
                      dam='graded',
                      n_i=50,
                      n_p=1000,
                      SIMPLIFY=TRUE)
stopCluster(cl)


results <- prep.results(control, sim_out)

saveRDS(results, 'sim_data/sim_2_graded_graded.rds')

results |>
  select(b, starts_with('omega')) |>
  pivot_longer(-b, names_to='IMV', values_to='value') |>
  group_by(b, IMV) |>
  summarize(mean = mean(value),
            se = sd(value),
            .groups='drop') |>
  ggplot(aes(x=b, y=mean, color=IMV,
             fill=IMV)) +
  geom_vline(aes(xintercept=-1), lty=2, alpha=0.5) +
  geom_vline(aes(xintercept=0), lty=2, alpha=0.5) +
  geom_vline(aes(xintercept=1), lty=2, alpha=0.5) +
  geom_line() +
  ggtitle('DGM: GRM - DAM: GRM') +
  theme_bw()

ggsave('fig/2_b_scan_graded_graded.png', width=8, height=6)

cl <- makeCluster(n_cores)
clusterSetRNGStream(cl, 1337)
clusterExport(cl, c('sim.fun', 'imv', 'imv_c', 'imv_t'))
sim_out <- clusterMap(cl, sim.fun,
                      b=control$b,
                      b_grid=c(-1,0,1),
                      dgm='graded',
                      dam='gpcm',
                      n_i=50,
                      n_p=1000,
                      SIMPLIFY=TRUE)
stopCluster(cl)


results <- prep.results(control, sim_out)

saveRDS(results, 'sim_data/sim_2_graded_gpcm.rds')

results |>
  select(b, starts_with('omega')) |>
  pivot_longer(-b, names_to='IMV', values_to='value') |>
  group_by(b, IMV) |>
  summarize(mean = mean(value),
            se = sd(value),
            .groups='drop') |>
  ggplot(aes(x=b, y=mean, color=IMV,
             fill=IMV)) +
  geom_vline(aes(xintercept=-1), lty=2, alpha=0.5) +
  geom_vline(aes(xintercept=0), lty=2, alpha=0.5) +
  geom_vline(aes(xintercept=1), lty=2, alpha=0.5) +
  geom_line() +
  ggtitle('DGM: GRM - DAM: GPCM') +
  theme_bw()

ggsave('fig/2_b_scan_graded_gpcm.png', width=8, height=6)


# Again fix threshold parameters at (-1, 0, 1) and vary two others along the
# interval [-2,2]. Make one plot per  as a function of both b parameters, with
# color indicating the value of omega. Include dashed horizontal and vertical lines
# at the values of the fixed threshold parameters.
