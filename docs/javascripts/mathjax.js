window.MathJax = {
  tex: {
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex"
  }
};

if (typeof document$ !== "undefined") {
  document$.subscribe(() => {
    if (
      typeof MathJax !== "undefined" &&
      MathJax.startup &&
      MathJax.startup.output
    ) {
      MathJax.startup.output.clearCache();
      MathJax.typesetClear();
      MathJax.texReset();
      MathJax.typesetPromise();
    }
  });
}
