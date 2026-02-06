import { StrictMode } from 'react'
import { createRoot } from 'react-dom/client'
import { CssBaseline, ThemeProvider, createTheme } from '@mui/material'
import './index.css'
import App from './App.tsx'

const theme = createTheme({
  palette: {
    mode: 'dark',
    background: {
      default: '#050816',
      paper: 'rgba(15, 23, 42, 0.92)',
    },
    primary: {
      main: '#38bdf8',
      light: '#7dd3fc',
      dark: '#0284c7',
      contrastText: '#000',
    },
    secondary: {
      main: '#a855f7',
      light: '#c084fc',
      dark: '#7c3aed',
      contrastText: '#000',
    },
    success: {
      main: '#22c55e',
      light: '#4ade80',
      dark: '#16a34a',
    },
    warning: {
      main: '#eab308',
      light: '#facc15',
      dark: '#ca8a04',
    },
    error: {
      main: '#ef4444',
      light: '#f87171',
      dark: '#dc2626',
    },
    text: {
      primary: '#e5e7eb',
      secondary: 'rgb(148, 163, 184)',
    },
    divider: 'rgba(148, 163, 184, 0.2)',
  },
  typography: {
    fontFamily: '"Inter", system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif',
    h2: {
      fontWeight: 700,
      letterSpacing: '-0.02em',
    },
    h4: {
      fontWeight: 700,
      letterSpacing: '0.5px',
    },
    h6: {
      fontWeight: 600,
      letterSpacing: '0.4px',
    },
    subtitle1: {
      fontWeight: 600,
      letterSpacing: '0.2px',
    },
    subtitle2: {
      fontWeight: 500,
      letterSpacing: '0.1px',
    },
    body2: {
      lineHeight: 1.6,
    },
    caption: {
      letterSpacing: '0.3px',
    },
  },
  shape: {
    borderRadius: 12,
  },
  components: {
    MuiCssBaseline: {
      styleOverrides: {
        body: {
          scrollbarColor: 'rgba(148, 163, 184, 0.3) rgba(15, 23, 42, 0.5)',
        },
      },
    },
    MuiButton: {
      styleOverrides: {
        root: {
          textTransform: 'none',
          fontWeight: 600,
        },
        contained: {
          boxShadow: '0 4px 14px rgba(56, 189, 248, 0.25)',
          '&:hover': {
            boxShadow: '0 6px 20px rgba(56, 189, 248, 0.35)',
          },
        },
      },
    },
    MuiPaper: {
      styleOverrides: {
        root: {
          backgroundImage: 'none',
        },
      },
    },
    MuiChip: {
      styleOverrides: {
        root: {
          fontWeight: 500,
        },
      },
    },
    MuiTextField: {
      styleOverrides: {
        root: {
          '& .MuiOutlinedInput-root': {
            '& fieldset': {
              borderColor: 'rgba(148, 163, 184, 0.3)',
            },
            '&:hover fieldset': {
              borderColor: 'rgba(148, 163, 184, 0.5)',
            },
            '&.Mui-focused fieldset': {
              borderColor: '#38bdf8',
            },
          },
        },
      },
    },
    MuiTooltip: {
      styleOverrides: {
        tooltip: {
          backgroundColor: 'rgba(15, 23, 42, 0.95)',
          border: '1px solid rgba(148, 163, 184, 0.2)',
          borderRadius: 8,
          fontSize: '0.75rem',
          padding: '8px 12px',
        },
        arrow: {
          color: 'rgba(15, 23, 42, 0.95)',
        },
      },
    },
    MuiTab: {
      styleOverrides: {
        root: {
          textTransform: 'none',
          fontWeight: 500,
          minHeight: 48,
        },
      },
    },
    MuiTabs: {
      styleOverrides: {
        indicator: {
          height: 3,
          borderRadius: '3px 3px 0 0',
        },
      },
    },
  },
})

createRoot(document.getElementById('root')!).render(
  <StrictMode>
    <ThemeProvider theme={theme}>
      <CssBaseline />
      <App />
    </ThemeProvider>
  </StrictMode>,
)
