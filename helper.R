# helper.R - Quarto 文獻驗證工具
# 功能：
#   1. 檢查 .qmd 檔案中的引用是否存在於 .bib 檔案
#   2. 驗證 .bib 檔案中的 DOI 是否有效

# 載入必要套件
if (!requireNamespace("httr", quietly = TRUE)) {
  install.packages("httr")
}
if (!requireNamespace("bib2df", quietly = TRUE)) {
  install.packages("bib2df")
}

library(httr)
library(bib2df)

#' 從 .qmd 檔案中提取所有引用 key
#' @param qmd_file .qmd 檔案路徑
#' @return 引用 key 的字元向量
extract_citations <- function(qmd_file) {
  if (!file.exists(qmd_file)) {
    stop("找不到檔案: ", qmd_file)
  }

  content <- readLines(qmd_file, warn = FALSE) |> paste(collapse = "\n")

  # 匹配 @citekey 格式（包含 [@key]、@key、[@key1; @key2] 等）
  pattern <- "@([a-zA-Z0-9_:.-]+)"
  matches <- gregexpr(pattern, content, perl = TRUE)
  keys <- regmatches(content, matches)[[1]]

  # 移除 @ 符號
  keys <- gsub("^@", "", keys)

  # 排除常見的非引用 @ 符號（如 email）
  keys <- keys[!grepl("\\.", keys) | grepl("^[a-zA-Z]+\\d{4}", keys)]

  unique(keys)
}

#' 從 .bib 檔案中提取所有 entry key
#' @param bib_file .bib 檔案路徑
#' @return 包含 key 和 doi 的 data.frame
parse_bib_file <- function(bib_file) {
  if (!file.exists(bib_file)) {
    stop("找不到檔案: ", bib_file)
  }

  tryCatch(
    {
      bib_df <- bib2df(bib_file)

      # 標準化欄位名稱（bib2df 會將欄位轉為大寫）
      names(bib_df) <- toupper(names(bib_df))

      result <- data.frame(
        key = bib_df$BIBTEXKEY,
        doi = if ("DOI" %in% names(bib_df)) bib_df$DOI else NA_character_,
        title = if ("TITLE" %in% names(bib_df)) bib_df$TITLE else NA_character_,
        stringsAsFactors = FALSE
      )

      return(result)
    },
    error = function(e) {
      # 備用方案：手動解析
      message("bib2df 解析失敗，使用備用解析器...")
      parse_bib_manual(bib_file)
    }
  )
}

#' 手動解析 .bib 檔案（備用方案）
#' @param bib_file .bib 檔案路徑
#' @return 包含 key 和 doi 的 data.frame
parse_bib_manual <- function(bib_file) {
  content <- readLines(bib_file, warn = FALSE) |> paste(collapse = "\n")

  # 提取 entry key
  key_pattern <- "@\\w+\\{([^,]+),"
  key_matches <- gregexpr(key_pattern, content, perl = TRUE)
  keys <- regmatches(content, key_matches)[[1]]
  keys <- gsub("@\\w+\\{|,", "", keys)
  keys <- trimws(keys)

  # 提取 DOI
  doi_pattern <- "doi\\s*=\\s*[{\"']?([^,\"'}+)[\"']?"

  # 分割成各個 entry
  entries <- strsplit(content, "(?=@\\w+\\{)", perl = TRUE)[[1]]
  entries <- entries[nchar(trimws(entries)) > 0]

  dois <- sapply(
    entries,
    function(entry) {
      match <- regexpr(doi_pattern, entry, ignore.case = TRUE, perl = TRUE)
      if (match > 0) {
        doi_text <- regmatches(entry, match)
        doi <- gsub(
          "doi\\s*=\\s*[{\"']?|[\"']?$",
          "",
          doi_text,
          ignore.case = TRUE
        )
        return(trimws(doi))
      }
      return(NA_character_)
    },
    USE.NAMES = FALSE
  )

  data.frame(
    key = keys,
    doi = dois[seq_along(keys)],
    title = NA_character_,
    stringsAsFactors = FALSE
  )
}

#' 驗證 DOI 是否有效
#' @param doi DOI 字串
#' @param timeout 請求超時秒數
#' @return 邏輯值，TRUE 表示有效
validate_doi <- function(doi, timeout = 10) {
  if (is.na(doi) || doi == "" || is.null(doi)) {
    return(NA)
  }

  # 清理 DOI
  doi <- trimws(doi)
  doi <- gsub("^https?://doi.org/", "", doi)
  doi <- gsub("^doi:", "", doi, ignore.case = TRUE)

  url <- paste0("https://doi.org/", doi)

  tryCatch(
    {
      response <- HEAD(url, timeout(timeout), config(followlocation = TRUE))
      return(status_code(response) == 200)
    },
    error = function(e) {
      return(FALSE)
    }
  )
}

#' 驗證所有文獻
#' @param qmd_files .qmd 檔案路徑（可多個）
#' @param bib_file .bib 檔案路徑
#' @param check_doi 是否驗證 DOI
#' @param verbose 是否顯示詳細資訊
#' @return 驗證結果的 list
validate_references <- function(
  qmd_files = NULL,
  bib_file = "references.bib",
  check_doi = TRUE,
  verbose = TRUE
) {
  # 自動尋找檔案
  if (is.null(qmd_files)) {
    qmd_files <- list.files(pattern = "\\.qmd$", recursive = TRUE)
    if (length(qmd_files) == 0) {
      stop("找不到任何 .qmd 檔案")
    }
  }

  if (!file.exists(bib_file)) {
    # 嘗試尋找 .bib 檔案
    bib_files <- list.files(pattern = "\\.bib$", recursive = TRUE)
    if (length(bib_files) > 0) {
      bib_file <- bib_files[1]
      if (verbose) message("使用找到的 .bib 檔案: ", bib_file)
    } else {
      stop("找不到 .bib 檔案")
    }
  }

  if (verbose) {
    message("=== Quarto 文獻驗證工具 ===\n")
    message("檢查的 .qmd 檔案: ", paste(qmd_files, collapse = ", "))
    message("使用的 .bib 檔案: ", bib_file, "\n")
  }

  # 提取所有引用
  all_citations <- character(0)
  for (qmd in qmd_files) {
    citations <- extract_citations(qmd)
    all_citations <- c(all_citations, citations)
    if (verbose) {
      message("從 ", qmd, " 提取到 ", length(citations), " 個引用")
    }
  }
  all_citations <- unique(all_citations)

  if (verbose) {
    message("\n總共 ", length(all_citations), " 個唯一引用\n")
  }

  # 解析 .bib 檔案
  bib_data <- parse_bib_file(bib_file)
  bib_keys <- bib_data$key

  if (verbose) {
    message(".bib 檔案中有 ", length(bib_keys), " 筆文獻\n")
  }

  # 檢查缺少的引用
  missing_citations <- setdiff(all_citations, bib_keys)
  unused_entries <- setdiff(bib_keys, all_citations)

  # 結果整理
  results <- list(
    total_citations = length(all_citations),
    total_bib_entries = length(bib_keys),
    missing_citations = missing_citations,
    unused_entries = unused_entries,
    doi_validation = NULL
  )

  # 輸出缺少的引用
  if (verbose) {
    message("--- 驗證結果 ---\n")

    if (length(missing_citations) > 0) {
      message("❌ 找不到的引用 (", length(missing_citations), " 個):")
      for (key in missing_citations) {
        message("   - @", key)
      }
      message("")
    } else {
      message("✅ 所有引用都存在於 .bib 檔案中\n")
    }

    if (length(unused_entries) > 0) {
      message("⚠️  未使用的文獻 (", length(unused_entries), " 個):")
      for (key in unused_entries) {
        message("   - ", key)
      }
      message("")
    }
  }

  # DOI 驗證
  if (check_doi) {
    if (verbose) {
      message("--- DOI 驗證 ---\n")
      message("正在驗證 DOI（這可能需要一些時間）...\n")
    }

    doi_results <- data.frame(
      key = bib_data$key,
      doi = bib_data$doi,
      valid = NA,
      stringsAsFactors = FALSE
    )

    has_doi <- !is.na(bib_data$doi) & bib_data$doi != ""

    if (verbose) {
      message("有 DOI 的文獻: ", sum(has_doi), " / ", nrow(bib_data), "\n")
    }

    # 驗證每個 DOI
    pb <- if (verbose && sum(has_doi) > 0) {
      txtProgressBar(min = 0, max = sum(has_doi), style = 3)
    } else {
      NULL
    }

    j <- 0
    for (i in which(has_doi)) {
      doi_results$valid[i] <- validate_doi(bib_data$doi[i])
      j <- j + 1
      if (!is.null(pb)) {
        setTxtProgressBar(pb, j)
      }
      Sys.sleep(0.5) # 避免請求過快
    }

    if (!is.null(pb)) {
      close(pb)
      message("")
    }

    results$doi_validation <- doi_results

    # 輸出 DOI 驗證結果
    if (verbose) {
      invalid_dois <- doi_results[
        !is.na(doi_results$valid) &
          doi_results$valid == FALSE,
      ]
      missing_dois <- doi_results[
        is.na(doi_results$doi) |
          doi_results$doi == "",
      ]

      if (nrow(invalid_dois) > 0) {
        message("\n❌ 無效的 DOI (", nrow(invalid_dois), " 個):")
        for (i in seq_len(nrow(invalid_dois))) {
          message("   - ", invalid_dois$key[i], ": ", invalid_dois$doi[i])
        }
      } else if (sum(has_doi) > 0) {
        message("\n✅ 所有 DOI 都有效")
      }

      if (nrow(missing_dois) > 0) {
        message("\n⚠️  缺少 DOI 的文獻 (", nrow(missing_dois), " 個):")
        for (key in missing_dois$key) {
          message("   - ", key)
        }
      }
    }
  }

  if (verbose) {
    message("\n=== 驗證完成 ===")
  }

  invisible(results)
}

#' 快速驗證（簡化版）
#' @param ... 傳遞給 validate_references 的參數
check_refs <- function(...) {
  validate_references(...)
}

# 如果直接執行此腳本
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) == 0) {
    # 使用預設值
    validate_references()
  } else {
    # 使用命令列參數
    qmd_files <- args[grepl("\\.qmd$", args)]
    bib_file <- args[grepl("\\.bib$", args)]

    if (length(qmd_files) == 0) {
      qmd_files <- NULL
    }
    if (length(bib_file) == 0) {
      bib_file <- "references.bib"
    }

    validate_references(qmd_files = qmd_files, bib_file = bib_file[1])
  }
}

# 使用說明
# ========
#
# 1. 基本使用（自動尋找檔案）：
#    source("helper.R")
#    check_refs()
#
# 2. 指定檔案：
#    validate_references(
#      qmd_files = c("paper.qmd", "appendix.qmd"),
#      bib_file = "references.bib"
#    )
#
# 3. 只檢查引用，不驗證 DOI：
#    validate_references(check_doi = FALSE)
#
# 4. 從命令列執行：
#    Rscript helper.R
#    Rscript helper.R paper.qmd references.bib served as a helper script to validate the citaiton
