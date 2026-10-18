// Renders the Mermaid diagrams of the README and of the function map article
// to PNG images, because the GitHub app does not draw Mermaid blocks. Run from
// the package root:
//   npm install --prefix <dir> mermaid@10.9.1
//   NODE_PATH=$(npm root -g) MERMAID_JS=<dir>/node_modules/mermaid/dist/mermaid.min.js \
//     node dev/diagrams/render.js
// Playwright and Chromium come from the global node installation. The widths
// printed are the natural widths in CSS pixels (the PNGs are drawn at 2x);
// the img tags of the README and the article use them.
const fs = require("fs");
const { chromium } = require("playwright");

const diagrams = [
  ["readme-function-map", "man/figures/README-function-map.png"],
  ["function-map-main", "vignettes/articles/function-map-main.png"],
  ["function-map-data", "vignettes/articles/function-map-data.png"],
  ["function-map-models", "vignettes/articles/function-map-models.png"],
  ["function-map-choice", "vignettes/articles/function-map-choice.png"],
  ["function-map-results", "vignettes/articles/function-map-results.png"],
];

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ deviceScaleFactor: 2, viewport: { width: 1400, height: 1400 } });
  // the label reset of the former pkgdown header, so that the boxes fit their text
  await page.setContent(`<html><head><style>
    body { margin: 0; background: #ffffff; font-family: "trebuchet ms", verdana, arial, sans-serif; }
    #d { display: inline-block; padding: 16px; background: #ffffff; }
    #d p { margin: 0; }
  </style></head><body><div id="d"></div></body></html>`);
  await page.addScriptTag({ path: process.env.MERMAID_JS });
  await page.evaluate(() => mermaid.initialize({ startOnLoad: false, securityLevel: "strict", theme: "default" }));
  for (const [name, out] of diagrams) {
    const src = fs.readFileSync(`dev/diagrams/${name}.mmd`, "utf8");
    await page.evaluate(async ([id, code]) => {
      const { svg } = await mermaid.render(id, code);
      const d = document.getElementById("d");
      d.innerHTML = svg;
      // draw the diagram at its natural size (mermaid scales it to the container)
      const el = d.querySelector("svg");
      const vb = el.viewBox.baseVal;
      el.setAttribute("width", vb.width);
      el.setAttribute("height", vb.height);
      el.style.maxWidth = "none";
    }, [name.replace(/-/g, "_"), src]);
    const box = await page.$("#d");
    await box.screenshot({ path: out });
    const size = await box.boundingBox();
    console.log(`${out}: width ${Math.round(size.width)} px`);
  }
  await browser.close();
})();
