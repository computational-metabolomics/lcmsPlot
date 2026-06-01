const path = require("path");

/** @type {import('tailwindcss').Config} */
module.exports = {
    content: [
        path.resolve(__dirname, "../../R/*.R")
    ],
    safelist: [
        // lp-* component classes are referenced from R templates.
        "lp-card",
        "lp-card-tight",
        "lp-section-title",
        "lp-btn",
        "lp-btn-primary",
        "lp-nav-item",
        "lp-nav-active",
        "lp-chip",
        "lp-divider",
        "lp-grid-main",
        "lp-grid-content",
        "lp-plot-card",
        // Defence-in-depth: regex-safelist every Tailwind utility we use
        // directly from R `class = "..."` strings. If the content scan
        // ever silently breaks again, these still end up in the output.
        { pattern: /^(flex|inline-flex|grid|block|hidden)$/ },
        { pattern: /^flex-(row|col)$/ },
        { pattern: /^items-(start|center|end)$/ },
        { pattern: /^justify-(start|center|end|between)$/ },
        { pattern: /^(gap|space-x|space-y)-(0|1|2|3|4|5|6|8)$/ },
        { pattern: /^(p|m|px|py|mx|my|mt|mb|ml|mr)-(0|1|2|3|4|5|6|8)$/ },
        { pattern: /^(w|h|min-h|min-w|max-h|max-w)-(full|screen|8|10|12)$/ },
        { pattern: /^rounded(-(sm|md|lg|xl|2xl|full))?$/ },
        { pattern:
            /^(bg|text|border)-(white|slate|indigo|brand)-(50|100|200|300|400|500|600|700|800)$/
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
                brand: {
                    50:  "#eef2ff",
                    100: "#e0e7ff",
                    500: "#6366f1",
                    600: "#4f46e5",
                    700: "#4338ca"
                }
            }
        }
    },
    plugins: []
};
