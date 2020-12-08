suppressPackageStartupMessages(require(BSgenome.Hsapiens.UCSC.hg38))
suppressPackageStartupMessages(require(furrr))
suppressPackageStartupMessages(require(future))
suppressPackageStartupMessages(require(glue))
suppressPackageStartupMessages(require(gpgr))
suppressPackageStartupMessages(require(here))
suppressPackageStartupMessages(require(kableExtra))
suppressPackageStartupMessages(require(tidyverse))

s <- tibble::tribble(
  ~sample, ~prefix, ~snv, ~sv, ~cnv, ~genome,
  "SBJ00037", "SBJ00037__SBJ00037_PRJ200529_L2000958", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00197", "SBJ00197__SBJ00197_MDX190186_L1900909", "-somatic-ensemble-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv", "GRCh37",
  "SBJ00481", "SBJ00481__SBJ00481_PTC_TsqN200327_L2000243", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00574_1", "SBJ00574_1__SBJ00574_PRJ200428_L2000802", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00574_2", "SBJ00574_2__SBJ00574_PRJ200429_L2000803", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00590",   "SBJ00590__SBJ00590_MDX200156_L2000853", "-somatic-PASS.vcf.gz",   "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00597",   "SBJ00597__SBJ00597_MDX200162_L2000855", "-somatic-PASS.vcf.gz",   "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00600",   "SBJ00600__SBJ00600_PRJ200484_L2000839", "-somatic-PASS.vcf.gz",   "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00602",   "SBJ00602__SBJ00602_PRJ200487_L2000841", "-somatic-PASS.vcf.gz",   "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00603",   "SBJ00603__SBJ00603_PRJ200489_L2000843", "-somatic-PASS.vcf.gz",   "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00605",   "SBJ00605__SBJ00605_PRJ200504_L2000847", "-somatic-PASS.vcf.gz",   "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00607",   "SBJ00607__SBJ00607_MDX200165_L2000857", "-somatic-PASS.vcf.gz",   "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00608", "SBJ00608__SBJ00608_MDX200170_L2000872", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00609", "SBJ00609__SBJ00609_PRJ200508_L2000874", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00610", "SBJ00610__SBJ00610_PRJ200511_L2000876", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00611", "SBJ00611__SBJ00611_PRJ200514_L2000878", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00612", "SBJ00612__SBJ00612_PRJ200517_L2000880", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00614", "SBJ00614__SBJ00614_PRJ200523_L2000916", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00615", "SBJ00615__SBJ00615_MDX200175_L2000910", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00616", "SBJ00616__SBJ00616_MDX200183_L2000912", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00617", "SBJ00617__SBJ00617_MDX200173_L2000908", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00618", "SBJ00618__SBJ00618_MDX200185_L2000914", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00620", "SBJ00620__SBJ00620_PRJ200532_L2000960", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00621", "SBJ00621__SBJ00621_MDX200196_L2000955", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00622", "SBJ00622__SBJ00622_MDX200190_L2000951", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00623", "SBJ00623__SBJ00623_MDX200193_L2000953", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00624", "SBJ00624__SBJ00624_MDX200199_L2000957", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00625", "SBJ00625__SBJ00625_PRJ200539_L2000962", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00627", "SBJ00627__SBJ00627_MDX200205_L2000998", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00628", "SBJ00628__SBJ00628_PRJ200545_L2000983", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00629", "SBJ00629__SBJ00629_PRJ200542_L2000981", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00630", "SBJ00630__SBJ00630_MDX200202_L2000979", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00631", "SBJ00631__SBJ00631_PRJ200548_L2000985", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00632", "SBJ00632__SBJ00632_PRJ200550_L2000987", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00633", "SBJ00633__SBJ00633_MDX200211_L2000991", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00635", "SBJ00635__SBJ00635_MDX200214_L2001000", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00636", "SBJ00636__SBJ00636_PRJ200552_L2001002", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00637", "SBJ00637__SBJ00637_PRJ200555_L2001004", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00638", "SBJ00638__SBJ00638_PRJ200576_L2001006", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
  "SBJ00639", "SBJ00639__SBJ00639_PRJ200579_L2001008", "-somatic-PASS.vcf.gz", "-manta.vcf.gz", ".purple.cnv.somatic.tsv", "hg38",
) %>%
  #dplyr::filter(sample %in% c("", "")) %>%
  dplyr::mutate(
    dd = file.path(here::here("nogit")),
    snv = file.path(dd, "snv", paste0(prefix, snv)),
    sv = file.path(dd, "sv", paste0(prefix, sv)),
    cnv = file.path(dd, "cnv", paste0(prefix, cnv))
  )

future::plan(multicore)
res <- seq_len(nrow(s)) %>%
  furrr::future_map(function(i) {
  #purrr::map(function(i) {
    cat(s$sample[i], "CHORD\n")
    chord_sv_df <- gpgr::chord_mantavcf2df(s$sv[i])
    chord <- gpgr::chord_run(
      vcf.snv = s$snv[i],
      df.sv = chord_sv_df,
      sv.caller = "manta",
      sample.name = s$sample[i],
      ref.genome = s$genome[i])

    cat(s$sample[i], "HRDetect\n")
    hrdetect <- gpgr::hrdetect_run(
      nm = s$sample[i],
      snvindel_vcf = s$snv[i],
      sv_vcf = s$sv[i],
      cnv_tsv = s$cnv[i],
      genome = s$genome[i],
      snvoutdir = file.path(here::here("nogit"), "results", "hrdetect", s$sample[i])
    )
    list(chord = chord$prediction,
         hrdetect = hrdetect)
  })

saveRDS(res, here::here("nogit/results/chord_hrdetect_2020-12-08.rds"))

res <- readRDS(here::here("nogit/results/chord_hrdetect_2020-12-08.rds"))

chord <- purrr::map(res, "chord") %>% dplyr::bind_rows()
hrdetect <- purrr::map(res, "hrdetect") %>% dplyr::bind_rows()

print(chord)
print(hrdetect)

s <- "SBJ00632"
chord %>%
  filter(sample == s) %>%
  t()

hrdetect %>%
  filter(sample == s) %>%
  t()
