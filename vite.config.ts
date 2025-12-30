
import { defineConfig, loadEnv } from 'vite';
import react from '@vitejs/plugin-react';

export default defineConfig(({ mode }) => {
  const env = loadEnv(mode, process.cwd(), '');
  return {
    plugins: [react()],
    define: {
      // Polyfill process.env so the Google GenAI SDK works in the browser
      'process.env': {
        API_KEY: JSON.stringify(env.API_KEY || env.VITE_API_KEY || '')
      }
    }
  };
});
