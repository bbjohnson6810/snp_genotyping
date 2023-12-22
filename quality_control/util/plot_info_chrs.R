library(ggplot2, quietly=T)
library(viridis, quietly=T)

# read in arguments: file locations
args <- commandArgs(TRUE)
info_all_file <- args[1]
info_0.9_file <- args[2]
out_prefix <- args[3]

# read in data from files
cat('Reading INFO scores for all SNPs... \n')
info_all <- read.table(info_all_file, sep='\t', header=T)

# change chromosome names (for easier visualization on plots)
info_all$CHROM[which(info_all$CHROM == 'chr1')] <- 'chr_01'
info_all$CHROM[which(info_all$CHROM == 'chr2')] <- 'chr_02'
info_all$CHROM[which(info_all$CHROM == 'chr3')] <- 'chr_03'
info_all$CHROM[which(info_all$CHROM == 'chr4')] <- 'chr_04'
info_all$CHROM[which(info_all$CHROM == 'chr5')] <- 'chr_05'
info_all$CHROM[which(info_all$CHROM == 'chr6')] <- 'chr_06'
info_all$CHROM[which(info_all$CHROM == 'chr7')] <- 'chr_07'
info_all$CHROM[which(info_all$CHROM == 'chr8')] <- 'chr_08'
info_all$CHROM[which(info_all$CHROM == 'chr9')] <- 'chr_09'
info_all$CHROM[which(info_all$CHROM == 'chr10')] <- 'chr_10'
info_all$CHROM[which(info_all$CHROM == 'chr11')] <- 'chr_11'
info_all$CHROM[which(info_all$CHROM == 'chr12')] <- 'chr_12'
info_all$CHROM[which(info_all$CHROM == 'chr13')] <- 'chr_13'
info_all$CHROM[which(info_all$CHROM == 'chr14')] <- 'chr_14'
info_all$CHROM[which(info_all$CHROM == 'chr15')] <- 'chr_15'
info_all$CHROM[which(info_all$CHROM == 'chr16')] <- 'chr_16'
info_all$CHROM[which(info_all$CHROM == 'chr17')] <- 'chr_17'
info_all$CHROM[which(info_all$CHROM == 'chr18')] <- 'chr_18'
info_all$CHROM[which(info_all$CHROM == 'chr19')] <- 'chr_19'
info_all$CHROM[which(info_all$CHROM == 'chr20')] <- 'chr_20'
info_all$CHROM[which(info_all$CHROM == 'chrX')] <- 'chr_X'
info_all$CHROM[which(info_all$CHROM == 'chrY')] <- 'chr_Y'
info_all$CHROM[which(info_all$CHROM == 'chrM')] <- 'chr_M'

# subset to keep high-quality INFO scores
info_clean <- info_all[info_all$INFO_SCORE > 0.9,]

# remove mitochondrial chromosomes - too small to plot
info_all <- info_all[info_all$CHROM != 'chr_M',]
info_clean <- info_clean[info_clean$CHROM != 'chr_M',]


# function to calculate mean INFO scores across sliding windows
# default window size: 1 MB
get_mean_info_scores <- function(info_df, window_size=1e5){

    # create an empty data frame to store the result
    window_df <- data.frame(chr = character(), start = integer(), end = integer(), info_mean = integer(), stringsAsFactors = FALSE)

    # iterate over unique chromosomes 
    for (chromosome in unique(info_df$CHROM)) {
        
      # subset info scores for the current chromosome
      chromosome_df <- info_df[info_df$CHROM == chromosome,]
      chromosome_snps <- chromosome_df$POS

      # calculate the start and end positions for each window
      start_positions <- seq(min(chromosome_snps), max(chromosome_snps), by = window_size)
      end_positions <- pmin(start_positions + window_size - 1, max(chromosome_snps))

      # calculate the mean info score within each window
      for (i in seq_along(start_positions)) {
        start <- start_positions[i]
        end <- end_positions[i]
        mean_info <- mean(chromosome_df[chromosome_snps >= start & chromosome_snps <= end,]$INFO_SCORE)

        # append the window information to the result data frame
        window_df <- rbind(window_df, data.frame(chr = chromosome, start = start, end = end, info_mean = mean_info))

      } # end of mutation loop
        
    } # end of chromosome loop
    
    # convert "chr" column to a factor with specified levels - to plot chromosomes in correct order
    chr_order <- c('chr_01','chr_02','chr_03','chr_04','chr_05','chr_06','chr_07','chr_08','chr_09',
               'chr_10','chr_11','chr_12','chr_13','chr_14','chr_15','chr_16','chr_17','chr_18',
                'chr_19','chr_20','chr_X','chr_Y')
    window_df$chr <- factor(window_df$chr, levels = rev(chr_order))

    return(window_df)

}


# function to plot info scores
plot_info_scores <- function(info_df, filtered=FALSE, window_size=1e5){

    # create ggplot heatmap with inverted y-axis
    if (filtered==FALSE){
            p1 <-
            ggplot(info_df, aes(x = start, xend = end, y = chr, yend = chr, color = info_mean)) +
            geom_segment(linewidth = 3.5) + theme_test() +
            # use complete 'mako' palette for unfiltered snps
            scale_colour_viridis_c(begin=0, end=1, option='mako') +
            ggtitle(paste0('Mean INFO scores: all SNPs (window size: ', window_size/1e6, 'MB)')) +
            xlab('Position (Mb)') +
            ylab('Chromosome') +
            theme(plot.title = element_text(size=14)) +
            theme(axis.title.y = element_text(margin = margin(r=20))) +
            theme(axis.title.x = element_text(margin = margin(t=20))) +
            theme(legend.position = c(0.84, 0.5)) +
            labs(color='Mean INFO Score') 
        }
        else{
            p1 <-
            ggplot(info_df, aes(x = start, xend = end, y = chr, yend = chr, color = info_mean)) +
            geom_segment(linewidth = 3.5) + theme_test() +
            # use top 10% of 'mako' palette for filtered snps (INFO > 0.9)
            scale_colour_viridis_c(begin=0.9, end=1, option='mako') +
            ggtitle(paste0('Mean INFO scores: filtered SNPs (> 0.9) (window size: ', window_size/1e6, 'MB)')) +
            xlab('Position (Mb)') +
            ylab('Chromosome') +
            theme(plot.title = element_text(size=14)) +
            theme(axis.title.y = element_text(margin = margin(r=20))) +
            theme(axis.title.x = element_text(margin = margin(t=20))) +
            theme(legend.position = c(0.84, 0.5)) +
            labs(color='Mean INFO Score') 
        }
                                       
    return(p1)

} # end of plotting funtion


# get mean info scores
cat('Calculating mean INFO scores... \n')
all_1e5 <- get_mean_info_scores(info_all, window_size=1e5)
clean_1e5 <- get_mean_info_scores(info_clean, window_size=1e5)
all_1e6 <- get_mean_info_scores(info_all, window_size=1e6)
clean_1e6 <- get_mean_info_scores(info_clean, window_size=1e6)

# save plots

# all INFO scores
cat('Saving plots for all SNPs \n')
png(paste0(out_prefix, '_INFO_all_map_0.1MB.png'),1000,800,res=130)
plot_info_scores(all_1e5, filtered=FALSE, window_size=1e5)
dev.off()
png(paste0(out_prefix, '_INFO_all_map_1MB.png'),1000,800,res=130)
plot_info_scores(all_1e6, filtered=FALSE, window_size=1e6)
dev.off()

# filtered INFO > 0.9
cat('Saving plots for high-quality SNPs \n')
png(paste0(out_prefix, '_INFO_0.9_map_0.1MB.png'),1000,800,res=130)
plot_info_scores(clean_1e5, filtered=TRUE, window_size=1e5)
dev.off()
png(paste0(out_prefix, '_INFO_0.9_map_1MB.png'),1000,800,res=130)
plot_info_scores(clean_1e6, filtered=TRUE, window_size=1e6)
dev.off()

