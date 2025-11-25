#' Synthetic Macroalgae Cover Field
#'
#' Generates a broad-scale synthetic field of macroalgae cover by combining
#' the projected hard and soft coral cover fields. Unlike corals, macroalgae
#' typically expands to occupy the remaining available space, so its cover
#' is calculated as the difference between total available area and coral cover.
#'
#' @title Synthetic Macroalgae Cover Field
#'
#' @param all_pts_effects_hcc A data.frame containing the projected hard coral cover.
#' @param all_pts_effects_sc A data.frame containing the projected soft coral cover.
#'
#' @return A list with:
#'   \itemize{
#'     \item `all_pts_effects_ma` – long-format data frame of projected macroalgae cover
#'   }
#'
#' @author Murray
#' @export
synthetic_field_ma <- function(all_pts_effects_hcc, all_pts_effects_sc) {
  ## Do all this on the link scale so that can use cumsum
  all_pts_effects_ma <- all_pts_effects_hcc |>
    dplyr::rename(HCC=Value) |> 
    dplyr::bind_cols(all_pts_effects_sc |>
                dplyr::select(SC=Value)) |>
    dplyr::mutate(Total_Avail = 0.8 - plogis(HCC) + plogis(SC),
      ## MA = Total_Avail*rbeta(n(), 2, 1),
      MA = Total_Avail,
      Value = qlogis(MA)) |>
    dplyr::select(-HCC, -SC, -Total_Avail, -MA)
  list(all_pts_effects_ma = all_pts_effects_ma)
}
