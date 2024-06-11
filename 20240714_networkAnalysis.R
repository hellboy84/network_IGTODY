install.packages("igraph")
library(igraph)

#データ読み込み
g5 <- read.csv("20240213net5.csv")

#ネットワークグラフ作成@windows環境
g <- graph_from_data_frame(g5, directed=TRUE)
plot(g, layout=layout_with_fr, vertex.label.cex=0.8, 
     vertex.frame.color="white", edge.arrow.width=0.6, 
     vertex.size=4.5, edge.arrow.size=0.1, edge.curved=0.15)

#macの場合 plot時にフォントを追加。plot外の「par(family = "HiraKakuProN-W3")」ではうまくいかない。
plot(g, layout=layout_with_fr, vertex.label.family="HiraKakuProN-W3",
     vertex.label.cex=0.8, vertex.frame.color="white", edge.arrow.width=0.6, 
     vertex.size=4.5, edge.arrow.size=0.1, edge.curved=0.15)

#デンドログラム作成
ceb <- cluster_edge_betweenness(g)
plot_dendrogram(ceb, mode="hclust")

#ネットワークグラフとデンドログラムの掛け合わせ
plot(ceb, g, layout=layout_with_fr, vertex.label.cex=0.8, 
     vertex.frame.color="white", edge.arrow.width=0.6, vertex.size=4.5, 
     edge.arrow.size=0.1, edge.curved=0.15)

#レイアウト固定。固定しないと毎回出力結果は変わる。
#デンドログラムの分割方法が違けどネットワークグラフは共通にしたい時などに必要
coord <- layout_with_fr(g)
plot(ceb, g, layout=coord, vertex.label.cex=0.8, 
     vertex.frame.color="white", edge.arrow.width=0.6, vertex.size=4.5, 
     edge.arrow.size=0.1, edge.curved=0.15)


#推移性(global)とクラスター係数(local)計算
transitivity(g, type = "global")
transitivity(g, type = "local")
transitivity(g, type = "average")

#次数カウント&ヒストグラム化&CSV出力
deg <- degree(g, mode="all")
deg_in <- degree(g, mode="in")
deg_out <- degree(g, mode="out")
hist(deg, breaks=14)
write.csv(deg, file="deg.csv")

#密度計算：どちらでも同じ結果
edge_density(g, loops=F)
#有向ネットワークの場合
ecount(g)/(vcount(g)*(vcount(g)-1)) 

#ネットワーク距離の計算とヒストグラム&CSV出力
network_distances <- distances(g)
hist(network_distances, breaks=10)
write.csv(network_distances, file="network_distances.csv")

#平均最短経路：有向性考慮
#無向グラフで計算すると，Fの結果と一致する。有向でTの方が値が下がる理由が不明……まぁFを採用する
mean_distance(g, directed=F)
mean_distance(g, directed=T)

#相互性
reciprocity(g, ignore.loops = TRUE)

#3点モチーフ数の計測
(count_motifs <- motifs(g, size=3))

#モチーフ一覧
#Macではwindowsではなくquartz()
#モチーフにカウントした数を書き込む
co <- matrix(c(1,1, 0,0, 2,0), ncol =2, byrow = TRUE)
windows()
par(mfrow = c(4,4), mar = c(0,0,0,0))
for (i in 1:16) {
  plot(graph_from_isomorphism_class(3, i-1, directed = TRUE),
       layout = co, vertex.color = "red", vertex.frame.color = NA,
       vertex.label = NA, edge.color = "darkgray",
       edge.arrow.size = 1.5, edge.width = 4, frame = TRUE, margin = 0.1)
  text(0,0, count_motifs[i], cex =4, col = "blue")}

#ヌルモデルとのモチーフ比較(Configulation model)
#元のネットワークのノード数と次数(出入)のパラメータを引き継ぐ1000個のランダムネットワークを作る
c_real <- motifs(g, size=3)
c_rand <- data.frame()
for(i in 1:1000){
  g_rand <- sample_degseq(deg_out, deg_in, m="simple")
  c_rand <- rbind(c_rand,motifs(g_rand,3))
}
names(c_rand)<-0:15
#元のネットワークがその平均からどれだけとおざかっているかのZ socreを出す
z_scores <- (c_real-colMeans(c_rand))/sqrt(colMeans(c_rand**2)-colMeans(c_rand)**2) 
z_scores

#他のネットワークグラフで同じことをした時に比較できるようにZ scoreをさらに標準化する。
z_scores_length <- sqrt(sum(z_scores^2, na.rm = TRUE))
z_scores_length
normalized_z_scores <- ifelse(is.na(z_scores), NA, z_scores / z_scores_length)
print(normalized_z_scores)
