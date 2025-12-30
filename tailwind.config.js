/** @type {import('tailwindcss').Config} */
export default {
  content: [
    "./index.html",
    "./**/*.{js,ts,jsx,tsx}",
  ],
  darkMode: 'class',
  theme: {
    extend: {
      colors: {
        lemos: {
          red: '#FF4B4B',
          dark: '#0E1117',
          card: '#262730',
          text: '#FAFAFA',
          sub: '#A0A0A0'
        }
      },
      animation: {
        progress: 'progress 1s ease-in-out infinite',
      },
      keyframes: {
        progress: {
          '0%': { transform: 'translateX(-100%)' },
          '100%': { transform: 'translateX(100%)' }
        }
      }
    },
  },
  plugins: [],
}