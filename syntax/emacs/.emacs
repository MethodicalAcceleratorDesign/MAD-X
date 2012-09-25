;; Enable syntax highlighting

(global-font-lock-mode t)
(setq font-lock-maximum-decoration t)


;; add madx highlighting
(add-to-list 'load-path "~/.emacs.d")
(autoload 'madx-mode "madx" "MADX-mode" t)
(setq auto-mode-alist (append '(("\\.madx$" . madx-mode))
  auto-mode-alist))
