
import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';

export default defineConfig({
  plugins: [react()],
  // Removed define: { 'process.env': {} } to allow process.env.API_KEY usage if environment is configured
});
