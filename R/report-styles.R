# The looks of a report. Each sets the CSS of the HTML report and the fonts,
# colours, and table fills of the PDF report; `rules` draws the tables with
# rules above and below them and under their header, as in print.
.report_styles <- list(
    plain = list(
        css = "
body { font-family: system-ui, -apple-system, 'Segoe UI', sans-serif;
  color: #222; background: #fff; margin: 0; }
main { max-width: 1100px; margin: 0 auto; padding: 16px; }
table { border-collapse: collapse; margin: 8px 0; }
th, td { padding: 4px 8px; border-bottom: 1px solid #ddd; text-align: left; }
th { background: #f4f4f4; }
.lp-fields td:first-child { font-weight: 600; }
.lp-report-page { border-top: 2px solid #ddd; margin-top: 32px; }
.lp-note { color: #666; font-style: italic; }
#lp-search { padding: 6px 8px; width: 100%; max-width: 320px; }
",
        font_family = "sans",
        text_colour = "black",
        heading_colour = "black",
        heading_face = "bold",
        title_face = "bold",
        header_colour = "black",
        header_face = "bold",
        header_fill = NA,
        row_fills = NA,
        rules = FALSE),
    clean = list(
        css = "
:root { --lp-text: #1f2328; --lp-muted: #656d76; --lp-border: #e3e6ea;
  --lp-subtle: #f6f8fa; --lp-accent: #2f5d8a; }
body { font-family: system-ui, -apple-system, 'Segoe UI', sans-serif;
  line-height: 1.5; color: var(--lp-text); background: #f7f8fa; margin: 0; }
main { max-width: 1120px; margin: 0 auto; padding: 32px 16px 64px; }
h1 { font-size: 1.7em; font-weight: 650; letter-spacing: -0.01em;
  margin: 0 0 20px; }
h2 { font-size: 1.2em; font-weight: 600; margin: 0 0 12px; }
h3 { font-size: 0.75em; font-weight: 600; text-transform: uppercase;
  letter-spacing: 0.06em; color: var(--lp-muted); margin: 16px 0 6px; }
#lp-index, .lp-report-page { background: #fff; border-radius: 10px;
  border: 1px solid var(--lp-border); padding: 20px 24px; margin: 0 0 24px; }
#lp-index > p { color: var(--lp-muted); margin-top: 0; }
table { border-collapse: collapse; margin: 4px 0;
  font-variant-numeric: tabular-nums; }
th, td { padding: 6px 12px; text-align: left;
  border-bottom: 1px solid var(--lp-border); }
th { font-size: 0.75em; font-weight: 600; text-transform: uppercase;
  letter-spacing: 0.04em; color: var(--lp-muted);
  background: var(--lp-subtle); }
tbody tr:nth-child(even) td { background: #fbfcfd; }
.lp-index-table tbody tr:hover td { background: #eef3f8; }
.lp-fields td:first-child { color: var(--lp-muted); }
a { color: var(--lp-accent); text-decoration: none; }
a:hover { text-decoration: underline; }
.lp-status { display: inline-block; padding: 1px 8px; border-radius: 999px;
  font-size: 0.8em; font-weight: 600; }
.lp-status-ok { background: #e6f4ea; color: #1e6b34; }
.lp-status-no-data { background: #eef0f2; color: #59636e; }
.lp-status-error { background: #fdecec; color: #a4282a; }
.lp-note { color: var(--lp-muted); font-style: italic; }
.lp-back { font-size: 0.9em; margin-bottom: 0; }
#lp-search { font: inherit; padding: 8px 12px; width: 100%; max-width: 320px;
  border: 1px solid var(--lp-border); border-radius: 6px; margin: 0 0 8px; }
#lp-search:focus { outline: 2px solid var(--lp-accent); outline-offset: 1px; }
",
        font_family = "sans",
        text_colour = "#1F2328",
        heading_colour = "#1F2328",
        heading_face = "bold",
        title_face = "bold",
        header_colour = "#656D76",
        header_face = "bold",
        header_fill = "#F1F3F5",
        row_fills = c("white", "#F8F9FA"),
        rules = FALSE),
    journal = list(
        css = "
body { font-family: Georgia, 'Times New Roman', serif; line-height: 1.55;
  color: #111; background: #fff; margin: 0;
  font-variant-numeric: lining-nums; }
main { max-width: 980px; margin: 0 auto; padding: 40px 16px 64px; }
h1 { font-weight: normal; font-size: 2em; margin: 0 0 24px; }
h2 { font-weight: normal; font-size: 1.35em; padding-bottom: 4px;
  border-bottom: 1px solid #111; }
h3 { font-weight: normal; font-variant: small-caps; letter-spacing: 0.05em;
  font-size: 1.05em; margin: 16px 0 4px; }
table { border-collapse: collapse; margin: 8px 0;
  border-top: 2px solid #111; border-bottom: 2px solid #111;
  font-variant-numeric: lining-nums tabular-nums; }
thead th { font-weight: normal; font-style: italic;
  border-bottom: 1px solid #111; }
th, td { padding: 4px 16px 4px 0; text-align: left; }
a { color: #111; text-decoration-thickness: 1px; }
.lp-report-page { margin-top: 56px; }
.lp-fields td:first-child { font-style: italic; }
.lp-status-no-data { font-style: italic; color: #555; }
.lp-status-error { font-style: italic; color: #8b1a1a; }
.lp-note { font-style: italic; color: #555; }
#lp-search { font: inherit; padding: 4px 0; width: 100%; max-width: 320px;
  border: none; border-bottom: 1px solid #111; background: transparent; }
",
        font_family = "serif",
        text_colour = "black",
        heading_colour = "black",
        heading_face = "plain",
        title_face = "plain",
        header_colour = "black",
        header_face = "italic",
        header_fill = NA,
        row_fills = NA,
        rules = TRUE)
)

# Layout and behaviour shared by every style.
.report_base_css <- function(font_size) {
    paste0("
body { font-size: ", font_size + 3, "px; }
img { max-width: 100%; height: auto; display: block; margin: 8px 0; }
.lp-index-table th { cursor: pointer; user-select: none; }
.lp-index-table th.lp-asc::after { content: ' \\25B2'; }
.lp-index-table th.lp-desc::after { content: ' \\25BC'; }
.lp-tables { display: flex; flex-wrap: wrap; column-gap: 32px; }
.lp-tables > div, .lp-scroll { overflow-x: auto; max-width: 100%; }
")
}
