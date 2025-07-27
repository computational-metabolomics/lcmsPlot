xcms_utils <- list(
  has_rt_alignment_been_performed = function(xcms_obj) {
    xcms::processHistory(xcms_obj) %>%
      lapply(function(p) p@type == "Retention time correction") %>%
      unlist() %>%
      any()
  },

  group_names = function(object, mzdec = 0, rtdec = 0) {
    feat_defs <- xcms::featureDefinitions(object)
    mzfmt <- paste("%.", mzdec, "f", sep = "")
    rtfmt <- paste("%.", rtdec, "f", sep = "")
    gnames <- paste("M", sprintf(mzfmt, feat_defs[,"mzmed"]), "T",
                    sprintf(rtfmt, feat_defs[,"rtmed"]), sep = "")

    if (any(dup <- duplicated(gnames))) {
      for (dupname in unique(gnames[dup])) {
        dupidx <- which(gnames == dupname)
        gnames[dupidx] <- paste(gnames[dupidx], seq(along = dupidx), sep = "_")
      }
    }

    return(gnames)
  },

  format_feature_identifiers = function(features, num_digits_rt = 0, num_digits_mz = 0) {
    features %>%
      # Extract decorations from 'name' (text after the first '_')
      mutate(idsDeco = stringr::str_extract(name, "_.*$")) %>%
      # Replace NA with empty string for idsDeco
      mutate(idsDeco = ifelse(is.na(idsDeco), "", idsDeco)) %>%
      # Create the custom name
      mutate(namecustom = make.unique(
        paste0("M", round(mz, num_digits_mz),
               "T", round(rt, num_digits_rt),
               idsDeco)
      )) %>%
      # Collapse peakidx column to a comma-separated string
      mutate(peakidx = sapply(peakidx, function(x) paste(x, collapse = ","))) %>%
      # Reorder columns
      relocate(name, namecustom) %>%
      select(-idsDeco)
  }
)
