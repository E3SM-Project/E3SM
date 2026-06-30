window.MathJax = {
  tex: {
    tags: 'all',
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true
  },
};

document$.subscribe(() => {
  MathJax.typesetPromise()
})

