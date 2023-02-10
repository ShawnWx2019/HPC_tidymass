#################################
#       Prj: Tidymass
#       Assignment: Pick peaking
#       Author: shawn
#       Date: Jun 1, 2022
#       Location: HENU
##################################
TEST = "FALSE"
options(stringsAsFactors = F)
options(warn = -1)
suppressMessages(if (!require('getopt')) BiocManager::install('getopt'))
suppressMessages(if (!require('crayon')) BiocManager::install('crayon'))

library(getopt)

msg_yes = green$bold$italic
msg_no = red$bold$italic
msg_run = blue$bold$italic$underline
msg_warning = yellow$bold$italic


# args --------------------------------------------------------------------

command=matrix(c(
  'help', 'h', 0, 'logic', 'help information',
  'POS', 'p', 1, 'character', 'Positive model .mzXML file path". ',
  'NEG', 'n', 1, 'character','Negative model .mzXML file path',
  'ppm', 'm', 1, 'double', 'Related to m/z\n                         fluctuation of m/z value (ppm) from scan to scan - depends on the mass spectrometer accuracy\n                         eg:10',
  'threads', 't', 1, 'integer', 'How many threads do you want to use.\n                          eg:20\n                          default:8',
  'show_output', 'o', 2, 'Logic', 'Need tic, bpc and rt_correction_plot or not\n                          default: TRUE',
  'p_min', 'a', 2, 'integer', 'Minimal of peakwidth for peak detection.\n                          eg:10\n                          default:10',
  'p_max', 'b', 2, 'integer', 'Max of peakwidth for peak detection.\n                          eg:10\n                          default:10',
  'min_fraction', 'f', 2, 'double', 'Related to Samples,\n                         to be valid, a group must be found in at least minFraction*n samples, with n=number of samples for each class of samples. A minFraction=0.5 corresponds to 50%.\n                         n=10, minFraction=0.5 => found in at least 5 samples\n                         default:0.5',
  'QC_tag', 'g', 2, 'character', 'File path of QCn.mzXML file placed.\n                         eg: QC.mzXML files placed in NEG/QC, set -g "QC", QC.mzXML files placed in NEG/qc, set -g "qc"',
  'snthresh', 's', 2, 'double', 'Related to intensity,\n                         signal/noise ratio threshold\n                         default: 5',
  'noise', 'N', 2, 'double', 'Related to intensity,\n                         each centroid must be greater than the “noise” value\n                         default:500'
),byrow = T, ncol = 5)
args = getopt(command)

## help information
if (!is.null(args$help)) {
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status=1)
}

## default value
if (is.null(args$POS)){
  message(msg_no("-p error \nPlease fill in the file path of POS files correctly!"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status = 1)
}
if (is.null(args$NEG)){
  message(msg_no("-n error \nPlease fill in the file path of NEG files correctly!"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status = 1)
}

if (is.null(args$QC_tag)){
  message(msg_no("-g error \nPlease fill in the file path of QC files correctly!"))
  message(msg_yes(paste(getopt(command, usage = T), "\n")))
  q(status = 1)
}

if (is.null(args$ppm)){
  args$ppm = 10
}

if (is.null(args$threads)){
  args$threads = 8
}

if (is.null(args$snthresh)){
  args$snthresh = 5
}

if (is.null(args$noise)){
  args$noise = 500
}

if (is.null(args$min_fraction)){
  args$min_fraction = 0.5
}


if (is.null(args$p_min)){
  args$p_min = 10
}

if (is.null(args$p_max)){
  args$p_max = 60
}

if (is.null(args$show_output)){
  args$show_output = TRUE
}


# test --------------------------------------------------------------------
library(tidymass)
##test
if (TEST == "TRUE") {
  T.POS = "../raw/NEG";
  T.NEG = "../raw/POS";
  T.ppm = 10;
  T.threads = 8;
  T.show_output = FALSE;
  T.p_min = 10;
  T.p_max = 60;
  T.min_fraction = 0.5;
  T.QC_tag = "QC";
  T.snthresh = 5;
  T.noise = 500
} else {
  T.POS = args$POS;
  T.NEG = args$NEG;
  T.ppm = args$ppm;
  T.threads = args$threads;
  T.show_output = args$show_output;
  T.p_min = args$p_min;
  T.p_max = args$p_max;
  T.min_fraction = args$min_fraction;
  T.QC_tag = args$QC_tag;
  T.snthresh = args$snthresh;
  T.noise = args$noise
}

# process data ------------------------------------------------------------

##> POS model
##>
message(msg_run("Step1. Peak peaking in positive model. This may take a long time. Please be patient!"))
process_data(
  path = T.POS,
  polarity = "positive",
  ppm = T.ppm,
  threads = T.threads,
  snthresh = T.snthresh,
  noise = T.noise,
  peakwidth = c(T.p_min,T.p_max),
  output_tic = T.show_output,
  output_bpc = T.show_output,
  output_rt_correction_plot = T.show_output,
  min_fraction = T.min_fraction,
  fill_peaks = FALSE,
  group_for_figure = T.QC_tag
)

##> NEG model
##>
message(msg_run("Step2. Peak peaking in negative model. This may take a long time. Please be patient!"))
process_data(
  path = T.NEG,
  polarity = "negative",
  ppm = T.ppm,
  threads = T.threads,
  snthresh = T.snthresh,
  noise = T.noise,
  peakwidth = c(T.p_min,T.p_max),
  output_tic = T.show_output,
  output_bpc = T.show_output,
  output_rt_correction_plot = T.show_output,
  min_fraction = T.min_fraction,
  fill_peaks = FALSE,
  group_for_figure = T.QC_tag
)

message(msg_yes("Peak picking steps fininsed. Please check result in /NEG/Result or /POS/Result"))

