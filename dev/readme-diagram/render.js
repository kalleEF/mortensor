// Renders the README diagram (function-map.mmd) to man/figures/README-function-map.png,
// because the GitHub mobile app does not draw Mermaid blocks. Run from the package root:
//   npm install --prefix <dir> mermaid@10.9.1   (the version the pkgdown site loads)
//   NODE_PATH=$(npm root -g) MERMAID_JS=<dir>/node_modules/mermaid/dist/mermaid.min.js \
//     node dev/readme-diagram/render.js
// Playwright and Chromium come from the global node installation.
const fs = require("fs");
const { chromium } = require("playwright");

(async () => {
  const src = fs.readFileSync("dev/readme-diagram/function-map.mmd", "utf8");
  const browser = await chromium.launch();
  const page = await browser.newPage({ deviceScaleFactor: 2, viewport: { width: 900, height: 1200 } });
  // the same label reset as in _pkgdown.yml, so that the boxes fit their text
  await page.setContent(`<html><head><style>
    body { margin: 0; background: #ffffff; font-family: "trebuchet ms", verdana, arial, sans-serif; }
    #d { display: inline-block; padding: 16px; background: #ffffff; }
    #d p { margin: 0; }
  </style></head><body><div id="d"></div></body></html>`);
  await page.addScriptTag({ path: process.env.MERMAID_JS });
  await page.evaluate(async (code) => {
    mermaid.initialize({ startOnLoad: false, securityLevel: "strict", theme: "default" });
    const { svg } = await mermaid.render("fm", code);
    const d = document.getElementById("d");
    d.innerHTML = svg;
    // draw the diagram at its natural size (mermaid scales it to the container)
    const el = d.querySelector("svg");
    const vb = el.viewBox.baseVal;
    el.setAttribute("width", vb.width);
    el.setAttribute("height", vb.height);
    el.style.maxWidth = "none";
  }, src);
  const box = await page.$("#d");
  await box.screenshot({ path: "man/figures/README-function-map.png" });
  const size = await box.boundingBox();
  console.log(`width ${Math.round(size.width)} px, height ${Math.round(size.height)} px (CSS pixels; PNG at 2x)`);
  await browser.close();
})();
