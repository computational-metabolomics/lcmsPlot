const path = require("path");

/** @type {import('tailwindcss').Config} */
module.exports = {
    content: [
        path.resolve(__dirname, "../../R/*.R")
    ],
    safelist: [
        // lp-* component classes are referenced from R templates.
        "lp-shell",
        "lp-sidebar",
        "lp-sidebar-group",
        "lp-sidebar-body",
        "lp-wordmark",
        "lp-main",
        "lp-topbar",
        "lp-page-title",
        "lp-page-subtitle",
        "lp-content",
        "lp-card",
        "lp-card-title",
        "lp-group-label",
        "lp-subheading",
        "lp-help",
        "lp-btn",
        "lp-btn-primary",
        "lp-link-btn",
        "lp-nav-item",
        "lp-nav-active",
        "lp-chip",
        "lp-chip-muted",
        "lp-divider",
        "lp-grid-content",
        "lp-plot-card",
        "lp-toolbar",
        "lp-segmented",
        "lp-sample-list",
        "lp-options",
        "is-busy",
        // Defence-in-depth: regex-safelist every Tailwind utility we use
        // directly from R `class = "..."` strings. If the content scan
        // ever silently breaks again, these still end up in the output.
        { pattern: /^(flex|inline-flex|grid|block|hidden)$/ },
        { pattern: /^flex-(row|col|wrap|1)$/ },
        { pattern: /^items-(start|center|end)$/ },
        { pattern: /^justify-(start|center|end|between)$/ },
        { pattern: /^(gap|space-x|space-y)-(0|1|2|3|4|5|6|8)$/ },
        { pattern: /^(p|m|px|py|mx|my|mt|mb|ml|mr)-(0|1|2|3|4|5|6|8)$/ },
        { pattern: /^(w|h|min-h|min-w|max-h|max-w)-(full|screen|8|10|12)$/ },
        { pattern: /^rounded(-(sm|md|lg|xl|2xl|full))?$/ },
        { pattern:
            /^(bg|text|border)-(white|gray|brand)-(25|50|100|200|300|400|500|600|700|800|900)$/
        },
        { pattern: /^text-(xs|sm|base|lg|xl|2xl)$/ },
        { pattern: /^font-(normal|medium|semibold|bold)$/ },
        { pattern: /^(border|border-(t|b|l|r))$/ },
        { pattern: /^shadow(-sm|-md)?$/ },
        { pattern: /^(sticky|top-0|z-10)$/ },
        { pattern: /^grid-cols-[1-4]$/ }
    ],
    theme: {
        extend: {
            colors: {
                // TailAdmin-style brand blue and neutral greys.
                brand: {
                    25:  "#f2f7ff",
                    50:  "#ecf3ff",
                    100: "#dde9ff",
                    200: "#c2d6ff",
                    300: "#9cb9ff",
                    400: "#7592ff",
                    500: "#465fff",
                    600: "#3641f5",
                    700: "#2a31d8"
                },
                gray: {
                    25:  "#fcfcfd",
                    50:  "#f9fafb",
                    100: "#f2f4f7",
                    200: "#e4e7ec",
                    300: "#d0d5dd",
                    400: "#98a2b3",
                    500: "#667085",
                    600: "#475467",
                    700: "#344054",
                    800: "#1d2939",
                    900: "#101828"
                }
            },
            fontFamily: {
                sans: ["Outfit", "ui-sans-serif", "system-ui", "-apple-system",
                       "Segoe UI", "Roboto", "Helvetica Neue", "Arial",
                       "sans-serif"]
            }
        }
    },
    plugins: []
};
