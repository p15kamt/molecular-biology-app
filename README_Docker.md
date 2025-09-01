# 🐳 Docker Deployment Guide
## Molecular Biology Analysis Application

### 📋 Prerequisites

- Docker Engine 20.10+
- Docker Compose 2.0+
- At least 8GB RAM available
- 4+ CPU cores recommended

### 🚀 Quick Start

#### Option 1: Docker Compose (Recommended)

```bash
# Clone the repository
git clone https://github.com/your-repo/molecular-biology-app.git
cd molecular-biology-app

# Build and run with Docker Compose
docker-compose up --build

# Access the application
# Open browser: http://localhost:8501
```

#### Option 2: Docker Build & Run

```bash
# Build the Docker image
docker build -t molecular-biology-app .

# Run the container
docker run -d \
  --name molecular_biology_app \
  -p 8501:8501 \
  -v $(pwd)/data:/app/data \
  -v $(pwd)/temp_files:/app/temp_files \
  --memory=8g \
  --cpus=4 \
  molecular-biology-app
```

### 🔧 Configuration

#### Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `STREAMLIT_SERVER_PORT` | 8501 | Application port |
| `STREAMLIT_SERVER_ADDRESS` | 0.0.0.0 | Server address |
| `PYTHONUNBUFFERED` | 1 | Python output buffering |

#### Volume Mounts

- `./data:/app/data` - Persistent data storage
- `./temp_files:/app/temp_files` - Temporary processing files

#### Resource Limits

- **Memory**: 8GB limit, 4GB reservation
- **CPU**: 4 cores allocated
- **Storage**: Depends on dataset size

### 🏥 Health Monitoring

The container includes health checks:

```bash
# Check container health
docker ps

# View health check logs
docker inspect --format='{{.State.Health.Status}}' molecular_biology_app
```

### 📊 Performance Tuning

#### For Large Datasets (>1M cells):

```yaml
# docker-compose.yml modifications
services:
  molecular-biology-app:
    mem_limit: 16g
    mem_reservation: 8g
    cpus: '8.0'
    environment:
      - STREAMLIT_SERVER_MAX_UPLOAD_SIZE=4096
```

#### For Production Deployment:

```yaml
# Add reverse proxy
services:
  nginx:
    image: nginx:alpine
    ports:
      - "80:80"
      - "443:443"
    volumes:
      - ./nginx.conf:/etc/nginx/nginx.conf:ro
    depends_on:
      - molecular-biology-app
```

### 🐛 Troubleshooting

#### Common Issues:

1. **Out of Memory Errors**
   ```bash
   # Increase memory limit
   docker-compose down
   # Edit docker-compose.yml: mem_limit: 16g
   docker-compose up
   ```

2. **Port Already in Use**
   ```bash
   # Change port mapping
   ports:
     - "8502:8501"  # Use different external port
   ```

3. **Permission Issues**
   ```bash
   # Fix volume permissions
   sudo chown -R $USER:$USER data temp_files
   chmod 755 data temp_files
   ```

### 📈 Monitoring & Logs

```bash
# View application logs
docker-compose logs -f molecular-biology-app

# Monitor resource usage
docker stats molecular_biology_app

# Access container shell
docker exec -it molecular_biology_app /bin/bash
```

### 🔄 Updates & Maintenance

```bash
# Update application
git pull origin main
docker-compose down
docker-compose up --build

# Clean up old images
docker image prune -f

# Backup data
tar -czf backup_$(date +%Y%m%d).tar.gz data/
```

### 🌐 Network Configuration

#### Custom Network:
```yaml
networks:
  molecular_biology_network:
    driver: bridge
    ipam:
      config:
        - subnet: 172.20.0.0/16
```

#### SSL/TLS Setup:
```bash
# Generate self-signed certificate
openssl req -x509 -newkey rsa:4096 -keyout key.pem -out cert.pem -days 365 -nodes

# Add to nginx configuration
server {
    listen 443 ssl;
    ssl_certificate /etc/ssl/cert.pem;
    ssl_certificate_key /etc/ssl/key.pem;
    
    location / {
        proxy_pass http://molecular-biology-app:8501;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
    }
}
```

### 📋 System Requirements Summary

| Component | Minimum | Recommended | Production |
|-----------|---------|-------------|------------|
| RAM | 4GB | 8GB | 16GB+ |
| CPU Cores | 2 | 4 | 8+ |
| Storage | 10GB | 50GB | 200GB+ |
| Docker | 20.10 | 24.0+ | Latest |

### 🆘 Support

For issues and questions:
1. Check logs: `docker-compose logs`
2. Review health status: `docker ps`
3. Verify resource usage: `docker stats`
4. Consult application logs in Streamlit interface

### 📚 Additional Resources

- [Docker Documentation](https://docs.docker.com/)
- [Docker Compose Reference](https://docs.docker.com/compose/)
- [Streamlit Deployment Guide](https://docs.streamlit.io/deploy)
- [Application Documentation](./ΤΕΧΝΙΚΟ_ΗΜΕΡΟΛΟΓΙΟ.md)
