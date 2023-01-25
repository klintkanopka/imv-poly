setwd('~/projects/imv-poly')
source('imvhelper.R')
source('simhelper.R')

library(tidyverse)
library(parallel)

n_cores <- 20

b_len <- 40
bs <- seq(from=-2, to=2, length.out=b_len)
N <- 100

control <- expand.grid(bs, bs) |>
  select(b1 = Var1, b2 = Var2) |>
  filter(b1 != b2)

control <- do.call(rbind, lapply(seq(N), function(k) cbind(id = k, control)))

cl <- makeCluster(n_cores)
clusterSetRNGStream(cl, 1337)
clusterExport(cl, c('sim.fun.2', 'imv', 'imv_c', 'imv_t'))
sim_out <- clusterMap(cl, sim.fun.2,
                      b1=control$b1,
                      b2=control$b2,
                      b_grid=c(-1,0,1),
                      dgm='graded',
                      dam='graded',
                      n_i=50,
                      n_p=1000,
                      SIMPLIFY=TRUE)
stopCluster(cl)


results <- prep.results(control, sim_out)

saveRDS(results, 'sim_data/sim_2a2_graded_graded.rds')

results |>
  group_by(b1, b2) |>
  summarize(across(starts_with('omega'), mean),
            .groups='drop') |>
  ggplot(aes(x=b1, y=b2, fill=omega_0)) +
  geom_tile() +
  geom_vline(aes(xintercept=-1), lty=2, alpha=0.5) +
  geom_vline(aes(xintercept=0), lty=2, alpha=0.5) +
  geom_vline(aes(xintercept=1), lty=2, alpha=0.5) +
  geom_hline(aes(yintercept=-1), lty=2, alpha=0.5) +
  geom_hline(aes(yintercept=0), lty=2, alpha=0.5) +
  geom_hline(aes(yintercept=1), lty=2, alpha=0.5) +
  ggtitle('DGM: GRM - DAM: GRM') +
  coord_equal() +
  theme_bw()

ggsave('fig/2_2b_scan_omega_0.png', width=8, height=6)

results |>
  group_by(b1, b2) |>
  summarize(across(starts_with('omega'), mean),
            .groups='drop') |>
  ggplot(aes(x=b1, y=b2, fill=omega_c)) +
  geom_tile() +
  geom_vline(aes(xintercept=-1), lty=2, alpha=0.5) +
  geom_vline(aes(xintercept=0), lty=2, alpha=0.5) +
  geom_vline(aes(xintercept=1), lty=2, alpha=0.5) +
  geom_hline(aes(yintercept=-1), lty=2, alpha=0.5) +
  geom_hline(aes(yintercept=0), lty=2, alpha=0.5) +
  geom_hline(aes(yintercept=1), lty=2, alpha=0.5) +
  ggtitle('DGM: GRM - DAM: GRM') +
  coord_equal() +
  theme_bw()

ggsave('fig/2_2b_scan_omega_c.png', width=8, height=6)

results |>
  group_by(b1, b2) |>
  summarize(across(starts_with('omega'), mean),
            .groups='drop') |>
  ggplot(aes(x=b1, y=b2, fill=omega_t)) +
  geom_tile() +
  geom_vline(aes(xintercept=-1), lty=2, alpha=0.5) +
  geom_vline(aes(xintercept=0), lty=2, alpha=0.5) +
  geom_vline(aes(xintercept=1), lty=2, alpha=0.5) +
  geom_hline(aes(yintercept=-1), lty=2, alpha=0.5) +
  geom_hline(aes(yintercept=0), lty=2, alpha=0.5) +
  geom_hline(aes(yintercept=1), lty=2, alpha=0.5) +
  ggtitle('DGM: GRM - DAM: GRM') +
  coord_equal() +
  theme_bw()

ggsave('fig/2_2b_scan_omega_t.png', width=8, height=6)
